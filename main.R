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
         time = as.POSIXlt(time, tz="America/New_York")) |>
  na.omit() |>
  as_tibble() 
  
         # skid = lead(real_open) - real_close,
skid <- 0.5
results <- tibble() 
colnames(results) <- unlist(str_split("j, slow_lag, fast_lag,
            ICAGR, drawdown, bliss, lake, end_val, trade_test", ", "))


######################## initiate optimization sequence
EMA_low <- 2
EMA_high <- 25
runs <- expand.grid(lag = seq(EMA_low, EMA_high, 1))
for (j in seq_len(nrow(runs))) {   
  df <- h4HA
  lag <- runs$lag[j]
  df$lag <- ewmaRcpp(df$real_close, runs$lag[j])
  
  df<- df |>                # create trade signal 
    mutate(cross = ha.close - lag,
           on = ifelse(cross > 0 & lag(cross) < 0, 1, 0),
           off = ifelse(cross < 0 & lag(cross) > 0, -1, 0))
  df <- slice(df, 30:n()) # drop first 30 rows to warm up EMA & Heikin Ashi
  
  if(df$cross[1] > 0) df$on[1] <- 1  # first row
  df <- df|>                         # trade signal details 
    mutate(signal = on + off,
          buy_date = if_else(on == 1, as_datetime(lead(time)), as.POSIXct(0)),
          buy_price = ifelse(on == 1, real_close + skid, 0),
          buy_amount = 0,
          sell_date = if_else(off == -1, as_datetime(lead(time)), as.POSIXct(0)),
          sell_price = ifelse(off == -1, real_close - skid, 0))

  if (df$cross[nrow(df)] > 0) {    # close out trade if long at EOF
  df$off[nrow(df)] <- -1
  df$signal[nrow(df)] <- -1
  df$sell_date[nrow(df)] <- df$time[nrow(df)]
  df$sell_price[nrow(df)] <- df$real_close[nrow(df)] - skid
}
 # return and risk variables
  start_value <- 1162 ; df$closed_pnl <- 1162 ; df$highwater <- 1162
  df$open_pnl <- 0 ; df$equity <- 0 ; df$water <- 0
  df$drawdown <- 0 ; buy_amount <- 0 ; buy_price <- 0
  # heat <- 0.1 # risk budget; not currently used, ATR_multiplier not used either
  
  # i=1 pnl trade calculation
  if(df$signal[1] == 1) {
    buy_amount = 1
    df$buy_amount[1] = buy_amount
    buy_price = df$buy_price[1]
  }
  # i>1 pnl trade calculation
  for (i in seq2(2, nrow(df))) {
    df$closed_pnl[i] = df$closed_pnl[i-1]
    if(df$signal[i] == -1) {
      df$closed_pnl[i] = df$closed_pnl[i] + 
        (df$sell_price[i] - buy_price) * buy_amount
      buy_amount = 0
      df$open_pnl[i] = 0
  
    } else if(df$signal[i] == 1) {
      buy_amount = 1
      df$buy_amount[i] = buy_amount
      buy_price = df$buy_price[i]
    } else {
      df$open_pnl[i] = (df$real_close[i] - buy_price) * buy_amount
    }
    df$equity[i] <- df$open_pnl[i] + df$closed_pnl[i]
    df$highwater[i] <- pmax(df$equity[i], df$highwater[i-1])
    df$water[i] <- df$highwater[i] - df$equity[i]
    df$drawdown[i] <- df$water[i] / df$highwater[i]
  }
  # summary trade table 
  trade_test <- sum(df$signal)
  if(trade_test == 0) {
    buys <- df |>
      select(on, buy_date, buy_price, buy_amount) |>
      filter(on == 1) |>
      select(!on)
    sells <- df |>
      select(off, sell_date:drawdown) |>
      filter(off == -1) |>
      select(!off)
    trades <- bind_cols(buys, sells)
    trades$buy_date <- as_datetime(as.numeric(trades$buy_date))
    trades$sell_date <- as_datetime(as.numeric(trades$sell_date))
    trades$trade_pnl <- (trades$sell_price - trades$buy_price) * trades$buy_amount
    # trades <- trades %>%
    #   mutate(min_low = map_dbl(1:n(), ~ {
    #     row <- .x
    #     df %>%
    #       filter(time > trades$buy_date[row], time <= trades$sell_date[row]) %>%
    #       pull(real_low) %>%
    #       min(na.rm = TRUE)
    #   }))
  }
end_value <- df$equity[nrow(df)] * 20
end_val <- end_value / 1000000
ratio <- end_value/ start_value
start_date <- min(df$time, na.rm=TRUE) ; end_date <- max(df$time, na.rm=TRUE)
date_range <- as.numeric(difftime(end_date, start_date, units = "days")) / 365.25
ICAGR <- if(ratio <= 0) 0 else log(ratio)/ date_range
drawdown <- max(df$drawdown)
lake <- sum(df$water) / sum(df$equity)
bliss <- ICAGR / drawdown
results[j,1:9] <- as_tibble_row(
  c(j=j, lag=lag, ICAGR=ICAGR, drawdown=drawdown, bliss=bliss, lake=lake, 
    end_value=end_value, end_val=end_val, trade_test=trade_test),
  .name_repair = "universal")

}       #################### optimization loop end
end_time <- Sys.time() ;forever <- end_time - start_time
secs <- forever  / nrow(runs)
sprintf("Yo, %1.2f total time and %1.4f per run, %i runs", forever, secs, nrow(runs))

# save the results
path <- getwd()
first_row_time <- df$time[1] ; second_row_time <- df$time[2]
interval_mins <- as.numeric(difftime(second_row_time, first_row_time, units = "mins"))
run_id <- paste0(get_month(start_date),"-", get_year(start_date),"--",
                 get_month(end_date), "-", get_year(end_date), " ",
                 interval_mins/60, "hr EMAs ", EMA_low, "-", EMA_high)
run_time <- paste0(" ", get_hour(end_time), ":", get_minute(end_time))
file_name <- paste0(path, "/results ", run_id, run_time, ".csv", sep="")
write_csv(results, file_name)

zz <- split_fun(trades, trade_pnl)
zz


################# the promised land of pretty graphs
options(ggrepel.max.overlaps = Inf)

df <- df |>
  mutate(date = as.POSIXct(time)) 
df |>        # The Market
  ggplot(aes(x = date, y = real_close)) +
  geom_line(size = 1, alpha = 0.5) + 
  labs(title=paste("The market"),
  subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr"))

df |>        # The Market, detailed  - maybe a different color?
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=real_low, ymax=real_high, x=date, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red") +
  geom_point(aes(x=date, y=real_close), shape=3, alpha=0.8) +
  labs(title=paste("The market: red is highs and lows, + are closes"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag ),) 

df |>        # Trades    Fix the damn trade lines and line color
  ggplot() +
  geom_point( aes(x=date, y=real_close), shape=3, alpha=0.4) +
  geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                    yend=sell_price, size = 2, color="black")) +
  labs(title="Trades and candle closes",
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag))

df |>        # lake over time with equity
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=equity*20, ymax=highwater*20, x=date, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red") +
  geom_line(aes(y = equity*20), size = 1, alpha = 0.8) +
  labs(title=paste("Lake Ratio with equity line"),
       subtitle=paste(start_date, "to", end_date, interval_mins/60, "hr", "EMA:", lag ),
       x="Year", y="Ending equity after $23k opening margin start") +
  scale_y_continuous(labels=scales::dollar_format())


df |>        # lake over time with highwter line
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=equity*20, ymax=highwater*20, x=date, fill = "band"), alpha = 0.9)+
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









