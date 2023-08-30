# short_trends trading identification and strategy calibration system
rm(list=ls())
library(tidyverse)
library(lubridate)
library(tidyquant)
library(quantmod)
library(Rcpp)
library(TTR)
library(plotly)
library(profvis)
library(microbenchmark)
library(rlang)
library(clock)
library(ggrepel)
library(zoo)

options(ggrepel.max.overlaps = Inf)      # ggrepel options for ggplot2
theme_set(theme_light())                # ggplot theme or _bw()

start_time <- Sys.time()

################ C code for EMA calculation with no leading NA
sourceCpp(
  code =
    "
     #include <Rcpp.h>
     // [[Rcpp::export]]
     Rcpp::NumericVector ewmaRcpp(Rcpp::NumericVector x, double a){
       int n = x.length();
       Rcpp::NumericVector s(n);
       s[0] = x[0];
       if (n > 1) {
         for (int i = 1; i < n; i++) {
           s[i] =  s[i-1] + (x[i] - s[i-1])/(( a + 1)/2);
         }
       }
       return s;
     }
    ")

#################### Heikin Ashi - 'average pace' also translates as 'average foot'
Rcpp::cppFunction('NumericMatrix RawHeikinAshi(NumericMatrix x, CharacterVector n) {
// assumes OHLC matrix input
int nrow = x.nrow(), ncol = 4, Op=0, Hi=1, Lo=2, Cl=3;
NumericMatrix ha(nrow,ncol);
for (int i = 0; i < nrow; i++) {
    ha(i, Cl) = (x(i,Op) + x(i,Hi) + x(i,Lo) + x(i,Cl)) / 4.0;
    ha(i, Op) = (i > 0) ? ((ha(i - 1, Op) + ha(i - 1, Cl)) / 2.0) : x(i, Op);
    ha(i, Hi) = std::max(x(i, Hi), std::max(ha(i, Op), ha(i, Cl)));
    ha(i, Lo) = std::min(x(i, Lo), std::min(ha(i, Op), ha(i, Cl)));
}
colnames(ha) = n;
return ha;
}')

HAOHLC <- function(x) {
  x <- OHLC(try.xts(x))
  r <- RawHeikinAshi(x, paste0("ha.", colnames(x)))
  return(reclass(r, x))
}

# calculates pnl summary statistics
split_fun <- function(data, column_name, factor_name) {
  summary_stats <- data |>
    mutate(category = ifelse({{ column_name }} >= 0, "Positive", "Negative")) |>
    group_by({{ factor_name }}, category) |>
    summarise(
      sum = sum({{ column_name }}),
      count = n(),
      mean = mean({{ column_name }}),
      median = median({{ column_name }}),
      sd = sd({{ column_name }}),
      min = min({{ column_name }}),
      max = max({{ column_name }})
    ) |>
    select({{ factor_name }}, category, sum, count, everything())
  return(summary_stats)
}

# Start dates: 1hr starts 1/1/20, 2hr 1/2/17, 3hr & 4hr 8/14/13, 6hr 1/2/04
# m30 <- read_csv("CME_MINI_NQ1!, 30_ac16d.csv", col_names = TRUE)
# # h1 <- read_csv("CME_MINI_NQ1!, 60_f371d.csv", col_names = TRUE)
# # h2 <- read_csv("CME_MINI_NQ1!, 120_925b9.csv", col_names = TRUE)
# # h3 <- read_csv("CME_MINI_NQ1!, 180_62ca4.csv", col_names = TRUE)
# # h6 <- read_csv("CME_MINI_NQ1!, 360_2a174.csv", col_names = TRUE)
# # spec(h6)
# NQ1D <- read_csv("NQ1D.csv", col_names = TRUE)
# NQ240 <- read_csv("NQ240.csv", col_names = TRUE)
h4 <- read_csv("Aug26b.csv", col_names = TRUE)

# h4 <- read_csv("CME_MINI_NQ1!, 240_48328.csv", col_names = TRUE)
HA_input <- select(h4, time:close)
h4HA <- HAOHLC(HA_input) 
h4HA <- h4HA |> rownames_to_column("time") |>
  mutate(real_open = h4$open,
         real_high = h4$high,
         real_low = h4$low,
         real_close = h4$close,
         time = as.POSIXct(time, tz="America/New_York"),
         yearmo = as.yearmon(time)) |>
  na.omit() |>
  as_tibble() 
  
first_row_time <- h4HA$time[1] ; second_row_time <- h4HA$time[2]
interval_mins <- as.numeric(difftime(second_row_time, first_row_time, units = "mins"))
results <- tibble() 
colnames(results) <- unlist(str_split("j, slow_lag, fast_lag,
            ICAGR, drawdown, bliss, lake, end_val, trade_test", ", "))

start_value <- 1162  # required overnight margin in points
skid <- 0.5   # skid is expected loss on trade execution
EMA_low <- 2
EMA_high <- 4
runs <- expand.grid(lag = seq(EMA_low, EMA_high, 1))
######################## initiate optimization sequence ########################
for (j in seq_len(nrow(runs))) {
  df <- h4HA
  lag <- runs$lag[j]
  df$lag <- ewmaRcpp(df$real_close, runs$lag[j])
  
  df<- df |>                # create trade signal 
    mutate(cross = ha.close - lag,
           on = if_else(cross > 0 & lag(cross) < 0, 1, 0),
           off = if_else(cross < 0 & lag(cross) > 0, -1, 0))
  
  df <- slice(df, 30:n()) # drop first 30 rows to warm up EMA & Heikin Ashi
  if(df$cross[1] > 0) df$on[1] <- 1  # first row only, set signal 
  
  df <- df|>                         # trade signal details 
    mutate(signal = on + off,
          buy_date = if_else(on == 1, as_datetime(lead(time), tz="America/New_York"), NA),
          sell_date = if_else(off == -1, as_datetime(lead(time), tz="America/New_York"), NA),
          buy_price = if_else(on == 1, real_close + skid, NA),
          sell_price = if_else(off == -1, real_close - skid, 0),
          buy_amount = 1,
          open_trade = cumsum(signal)) |>
    fill(buy_price) |>
    fill(buy_date) |>
    drop_na(buy_price)
  
  if (df$cross[nrow(df)] > 0 | df$signal[nrow(df)] == -1) { # close out trade if long at EOF
  df$off[nrow(df)] <- -1
  df$signal[nrow(df)] <- -1
  df$open_trade[nrow(df)] <- 0
  df$sell_date[nrow(df)] <- df$time[nrow(df)]
  df$sell_price[nrow(df)] <- df$real_close[nrow(df)] - skid       
}
  df <- df |>
    mutate(
      trade_pnl = if_else(signal == -1, (sell_price - buy_price) * buy_amount, 0),
      closed_pnl = cumsum(trade_pnl) + start_value,
      open_pnl = (real_close - buy_price) * buy_amount * open_trade,
      equity = open_pnl + closed_pnl,
      highwater = cummax(equity),
      lake = highwater - equity,
      drawdown = lake / highwater
    )

  # summary trade table 
  trade_test <- sum(df$signal) 
  if(trade_test != 0) warning("HOLY SHIT! THE TRADES ARE FUCKED!")
    trades <- df |>
      select(off, buy_date:drawdown) |>
      filter(off == -1) |>
      mutate(
        off=NULL, open_trade=NULL, buy_amount=NULL,
        color = ifelse(trade_pnl>0, "green", "red")) 

    ##################### make the MAE column   ######################################
        
    # trades <- trades %>%
    #   mutate(min_low = map_dbl(1:n(), ~ {
    #     row <- .x
    #     df %>%
    #       filter(time > trades$buy_date[row], time <= trades$sell_date[row]) %>%
    #       pull(real_low) %>%
    #       min(na.rm = TRUE)
    #   }))
    
end_value <- df$equity[nrow(df)] 
end_val <- end_value / 1000000
ratio <- end_value/ start_value
start_date <- min(df$time, na.rm=TRUE) 
end_date <- max(df$time, na.rm=TRUE)
date_range <- as.numeric(difftime(end_date, start_date, units = "days")) / 365.25
ICAGR <- if(ratio <= 0) 0 else log(ratio)/ date_range
drawdown <- max(df$drawdown)
lake <- sum(df$lake) / sum(df$equity)
bliss <- ICAGR / drawdown
results[j,1:9] <- as_tibble_row(          
  c(j=j, lag=lag, ICAGR=ICAGR, drawdown=drawdown, bliss=bliss, lake=lake, 
    end_value=end_value, end_val=end_val, trade_test=trade_test),
  .name_repair = "universal")

# save the results
path <- getwd()
run_id <- paste0(get_month(start_date),"-", get_year(start_date),"--",
                 get_month(end_date), "-", get_year(end_date), " ",
                 interval_mins/60, "hr EMAs ", EMA_low, "-", EMA_high)
run_time <- paste0(" ", get_hour(start_time), ":", get_minute(start_time))
file_name <- paste0(path, "/results ", run_id, run_time, ".csv", sep="")
write_csv(results, file_name)

df |>        # Trades, EMA and market graph  
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=real_low, ymax=real_high, x=time, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="lightblue") +
  geom_point(aes(x=time, y=real_close), shape=3, alpha=0.8) +
  geom_line(aes(x=time, y=lag), alpha=0.9) +
  geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                                yend=sell_price), linewidth = 2, color=trades$color) +
  labs(title=sprintf("NQ: %.0f trades, EMA: %0.f, ICAGR:, %.2f, bliss: %.2f, lake: %.2f", 
                     nrow(trades), lag, ICAGR, bliss, lake),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr" )) +
  xlab("Date")+
  ylab("NQ") 

ggsave(paste("output/fig 1b EMA", lag, interval_mins/60, "hr", as.Date(start_date), "to", as.Date(end_date), ".pdf"))



}       #################### optimization loop end    ##########################

end_time <- Sys.time() ;forever <- end_time - start_time
secs <- forever  / nrow(runs)
sprintf("Yo, %1.2f total time and %1.2f per run, %i runs", forever, secs, nrow(runs))


zz <- split_fun(trades, trade_pnl)
zz


################# the promised land of pretty graphs

# df |>        # The Market, EMA and trades
#   ggplot(aes(x = date, y = real_close)) +
#   geom_line(size = 1, alpha = 1) + 
#   geom_line(aes(x=time, y=lag), alpha=0.5) +
#   labs(title=paste("NQ futures closes,", lag, "period EMA, for", nrow(trades),"trades"),
#   subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr")) +
#   theme_light()
# ggsave("output/fig1.pdf")

df |>        # Trades and market graph  
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=real_low, ymax=real_high, x=time, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="lightblue") +
  geom_point(aes(x=time, y=real_close), shape=3, alpha=0.8) +
  geom_line(aes(x=time, y=lag), alpha=0.5) +
  geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                                yend=sell_price), linewidth = 2, color=trades$color) +
  labs(title=paste("NQ: blue range of highs and lows, + are closes,", nrow(trades), "trades"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag )) +
  xlab("Date")+
  ylab("NQ")
ggsave(paste("output/fig 1b EMA", lag, interval_mins/60, "hr", as.Date(start_date), "to", as.Date(end_date), ".pdf"))

df |>        # Trades    
  ggplot() +
  geom_line( aes(x=time, y=real_close), size=1, alpha=1) +
  geom_line(aes(x=time, y=lag), alpha=0.4) +
  geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                    yend=sell_price), linewidth = 2, color=trades$color) +
  labs(title="Trades, EMA and candle closes",
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

tradez <- trades |>        # Trades facet wrap attempt   Fix the damn trade lines and line color
  mutate(
    buy_date_day_before = buy_date - days(1),
    sell_date_day_after = sell_date + days(1),
    trade_index = row_number()  # Create an index variable for trade facets
  ) |>
  ggplot() +
  geom_line(data=df, aes(x=time, y=lag), alpha=0.6) +
  geom_line(data=df, aes(x=time, y=real_close), alpha=0.9) +
  geom_segment(data=tradez, aes(x=buy_date_day_before, y=buy_price, 
        xend=sell_date_day_after, yend=sell_price, size = 1, color="black")) +
  facet_wrap(~ trade_index, scales = "free_x") +
  labs(title="Trades, EMA and candle closes",
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))



df |>        # lake over time with equity
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=equity*20, ymax=highwater*20, x=time, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red") +
  geom_line(aes(y = equity*20), size = 1, alpha = 0.8) +
  labs(title=paste("Lake Ratio with equity line"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag ),
       x="Year", y="Ending equity after $23k opening margin start") +
  scale_y_continuous(labels=scales::dollar_format())


df |>        # lake over time with highwter line
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=equity*20, ymax=highwater*20, x=time, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red") +
  geom_line(aes(y = highwater*20), size = 1, alpha = 0.8) +
  labs(title=paste("Lake Ratio with highwater line"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag ),
       x="Year", y="Ending equity after $23k opening margin start") +
  scale_y_continuous(labels=scales::dollar_format())

df |>       # equity, highwater and market lines
  mutate(market = (real_close-min(real_close))* (max(highwater) / (max(real_close)-min(real_close)))) |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = equity),size = 1,  alpha = 0.9) +
  geom_line(aes(y=highwater), size=1,  alpha=0.4, color="black") +
  geom_line(aes(y=market), size=2, alpha=0.1) +
  labs(title=paste("Equity with highwater and the market"),
    subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

#   risk and return optimization scatterplots
rzlt <- results
 
rzlt |>         # labels for EMA numbers, little white boxes
  ggplot(aes(x = ICAGR, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Growth rate vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

rzlt |>
  ggplot(aes(x = lake, y = bliss, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
            min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Lake ratio vs bliss"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

rzlt |>
  ggplot(aes(x = lake, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Lake ratio vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

rzlt |>
  ggplot(aes(x = end_value, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  scale_x_continuous(labels=scales::dollar_format()) +
  labs(title=paste("Ending value vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag)) 

  # geom_smooth(method = "lm")

rzlt |>
  ggplot(aes(x = ICAGR, y = lake, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.1, "lines"), label.size=0.05, 
       min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Growth rate vs lake ratio"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag)) 

rzlt |>
  ggplot(aes(x = bliss, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
          min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Bliss vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag)) 









