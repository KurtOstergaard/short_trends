# short_trends EMA trading identification and strategy calibration system
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
start_time <- Sys.time()               # for internal monitoring
run_time <- paste0(" ", get_hour(start_time), "-", get_minute(start_time))
path <- getwd()

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
#################### Heikin Ashi - Japanese for 'average pace' also translates as average foot or leg
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
h4 <- read_csv("CME_MINI_NQ1!, 240_48328.csv", col_names = TRUE, show_col_types = FALSE)
# # h6 <- read_csv("CME_MINI_NQ1!, 360_2a174.csv", col_names = TRUE)
# # spec(h6)
# NQ1D <- read_csv("NQ1D.csv", col_names = TRUE)
# NQ240 <- read_csv("NQ240.csv", col_names = TRUE)
# h4 <- read_csv("Aug26b.csv", col_names = TRUE, show_col_types = FALSE)

HA_input <- select(h4, time:close)     # create the Haikin Ashi candlesticks
h4HA <- HAOHLC(HA_input) 
h4HA <- h4HA |> rownames_to_column("time") |>
  mutate(real_open = h4$open,          # add back the regular candlesticks
         real_high = h4$high,
         real_low = h4$low,
         real_close = h4$close,
         time = as.POSIXct(time, tz="America/New_York"),
         yearmo = as.yearmon(time)) |>     # this call was from the zoo package
  na.omit() |>
  as_tibble() 
  
# discern time interval from input file
first_row_time <- h4HA$time[1] ; second_row_time <- h4HA$time[2] ; third_row_time <- h4HA$time[3]
interval <- min(as.numeric(difftime(second_row_time, first_row_time, units = "mins")),
                as.numeric(difftime(third_row_time, second_row_time, units = "mins")))
interval <- if(floor(interval/60) == interval/60) interval/60 else interval

start_date <- min(h4HA$time, na.rm=TRUE) 
end_date <- max(h4HA$time, na.rm=TRUE)
date_range <- as.numeric(difftime(end_date, start_date, units = "days")) / 365.25

trades_global = tibble() # hey, why not create a global trades file by EMA?
results <- tibble() # create the file to collect the results of each run
# colnames(results) <- unlist(str_split("j, lag,
#             ICAGR, drawdown, bliss, lake, end_val, trade_test", ", "))
epoch <- paste0(get_month(start_date),"-", get_day(start_date), "-",
                get_year(start_date)-2000,"to", get_month(end_date), 
                "-", get_day(end_date), "-", get_year(end_date)-2000)

start_value <- 1162  # required overnight margin in points
skid <- 0.5   # skid is expected loss on trade execution
EMA_low <- 2
EMA_high <- 25
runs <- expand.grid(lag = seq(EMA_low, EMA_high, 1))

######################## EMA optimization sequence ########################

for (j in seq_len(nrow(runs))) {
  
  df <- h4HA
  lag <- runs$lag[j]
  df$lag <- ewmaRcpp(df$real_close, runs$lag[j])
  
  df<- df |>                # create trade signal 
    mutate(cross = ha.close - lag,
           on = if_else(cross > 0 & lag(cross) < 0, 1, 0),
           off = if_else(cross < 0 & lag(cross) > 0, -1, 0))
  
  df <- slice(df, 30:n()) # drop first 30 rows to warm up EMA & Heikin Ashi
  
  df$on[1] <-  if_else(df$cross[1] > 0, 1, 0)  # catch first row signal, if there
  
  df <- df|>                         # trade details 
    mutate(signal = on + off,
           buy_date = if_else(on == 1, as_datetime(lead(time), tz="America/New_York"), NA),
           sell_date = if_else(off == -1, as_datetime(lead(time), tz="America/New_York"), NA),
           buy_price = if_else(on == 1, real_close + skid, NA),
           sell_price = if_else(off == -1, real_close - skid, 0),
           buy_amount = 1,             # no work on trade sizing yet, important though it is
           open_trade = cumsum(signal)) |>
    fill(buy_price) |>
    fill(buy_date)  |>
    drop_na(buy_price)
  
  if(df$on[nrow(df)] == 1) {     # no new trades in last period
    df$on[nrow(df)] = 0 ; df$signal[nrow(df)] = 0
  } else if (df$cross[nrow(df)] >0 | df$signal[nrow(df)] == -1) { # close out trade if long at EOF
    df$off[nrow(df)] <- -1
    df$signal[nrow(df)] <- -1
    df$open_trade[nrow(df)] <- 0
    df$sell_date[nrow(df)] <- df$time[nrow(df)]
    df$sell_price[nrow(df)] <- df$real_close[nrow(df)] - skid       
  }
  df <- df |>
    mutate(
      trade_pnl = if_else(signal == -1, (sell_price - buy_price) * buy_amount, 0),
      win_lose = as.factor(if_else(trade_pnl>0, "winner", "loser")),
      closed_pnl = cumsum(trade_pnl) + start_value,
      open_pnl = (real_close - buy_price) * buy_amount * open_trade,
      equity = open_pnl + closed_pnl,
      highwater = cummax(equity),
      lake = highwater - equity,
      drawdown = lake / highwater)
  
  # summary trade table 
  trade_test <- sum(df$signal) 
  if(trade_test != 0) warning("HOLY SHIT! THE TRADES ARE FUCKED!")
  trades <- df |>
    select(off, buy_date:drawdown) |>
    filter(off == -1) |>
    mutate(
      off=NULL, open_trade=NULL, buy_amount=NULL)
  
  trades <- trades |>
    mutate(MAE = map_dbl(1:n(), ~ {
      row <- .x
      df |>
        filter(time >= trades$buy_date[row], time <= trades$sell_date[row]) |>
        pull(real_low) |>
        min(na.rm = TRUE)
    }), .after=win_lose,
    MAE_mark = MAE - buy_price,
    MAE_percent = (MAE / buy_price) - 1,
    trade_pnl_percent = trade_pnl / buy_price)
  
  # risk and return calculation
  end_value <- df$equity[nrow(df)] 
  end_val <- end_value / 1000000
  ratio <- end_value/ start_value
  ICAGR <- if(ratio <= 0) 0 else log(ratio)/ date_range
  drawdown <- max(df$drawdown)
  lake <- sum(df$lake) / sum(df$equity)
  bliss <- ICAGR / drawdown
  trade_count <- nrow(trades)
  trade_total_pnl <- sum(trades$trade_pnl)
  zz <- split_fun(trades, trade_pnl)
  wins <- zz[2,3] ; losses <- zz[1,3] ; won <- zz[2,2]; lost <- zz[1,2]
  win_rate <- wins / trade_count ; dollar_won <- -zz[2,4]/zz[1,4]
  results[j,1:17] <- as_tibble_row(          
    c(j=j, lag=lag, ICAGR=ICAGR, drawdown=drawdown, bliss=bliss, lake=lake, 
      end_value=end_value, end_val=end_val, trade_test=trade_test, 
      trade_count=trade_count, wins=wins, losses=losses, win_rate=win_rate,
      trade_total_pnl=trade_total_pnl, won=won, lost=lost, dollar_won=dollar_won),
    .name_repair = "universal")
  
  # accumulate the trades into the global trade table
  trade_tmp <- trades |>
    mutate(EMA = lag)
  trades_global <- trades_global |>
    bind_rows(trade_tmp)
  
  df |>        # Killer graph: trades by EMA with equity and market
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin=real_low, ymax=real_high, x=time, fill = "band"), alpha = 0.9)+
    scale_fill_manual("", values="gray80") +
    geom_point(aes(x=time, y=real_close), shape=3, alpha=0.8) +
    geom_line(aes(x=time, y=lag), alpha=0.9) +
    geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                        yend=sell_price, color=factor(win_lose)), linewidth = 2) +
    scale_color_manual(values= c("red", "green3")) +
    scale_y_continuous(sec.axis=sec_axis(~ . *(max(df$real_high)/min(df$real_low))-(min(df$real_low)) )) +
    geom_line(aes(x=time, y=equity *(max(real_high)/min(real_low)) + min(real_low)-min(equity)), size=2, alpha=0.7, color="deepskyblue") +
    labs(title=sprintf("NQ: EMA: %0.f, %.0f trades, ICAGR:, %.2f, bliss: %.2f, lake: %.2f", 
                       lag, nrow(trades), ICAGR, bliss, lake),
         subtitle=paste(interval, "hr", start_date, "to", end_date)) +
    xlab("Date")+
    ylab("NQ") +
    theme(legend.position = "none")
          
    ggsave(paste0("output/run ", interval, "hr EMA ", lag, " ", epoch, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)
          
trades |>
  ggplot(aes(MAE_percent, trade_pnl_percent, size=3, color=factor(win_lose))) +
  geom_point(shape=16) +
  scale_color_manual(values= c("red","green3")) + 
  labs(title=sprintf("NQ: EMA: %0.f, %.0f trades, ICAGR:, %.2f, bliss: %.2f, lake: %.2f", 
                     lag, nrow(trades), ICAGR, bliss, lake),
       subtitle=paste(interval, "hr", start_date, "to", end_date)) +
  xlab("Maximum Adverse Excursion") + 
  ylab("Trade P&L") +
  theme(legend.position = "none")

ggsave(paste0("output/stop ", interval, "hr EMA ", lag, " ", epoch, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

df |>        # lake over time with equity line
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin=equity*20, ymax=highwater*20, x=time, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red") +
  geom_line(aes(y = equity*20), size = 1, alpha = 0.8) +
  labs(title=paste("Lake Ratio with equity line"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag ),
       x="Year", y="Ending equity after $23k opening margin start") +
  scale_y_continuous(labels=scales::dollar_format())
ggsave(paste0("output/lake over time ", interval, "hr EMA ", lag, " ", epoch, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

}       #################### optimization loop end    ##########################

# save the results and trades_global files
run_id <- paste0( " ", interval, "hr EMAs ", EMA_low, "-", EMA_high, " from", epoch)
results_file_name <- paste0(path, "/output/results ", run_id, run_time, ".csv", sep="")
write_csv(results, results_file_name)
trade_file_name <- paste0(path, "/output/trades ", run_id, run_time, ".csv", sep="")
write_csv(trades_global, trade_file_name)


end_time <- Sys.time() ;forever <- end_time - start_time
secs <- forever  / nrow(runs)
sprintf("Yo, %1.2f total time and %1.2f per run, %i runs, over %1.2f years of data", 
        forever, secs, nrow(runs), date_range)

################# the promised land of pretty graphs


####################   risk and return optimization scatterplots

results |>         # labels for EMA numbers, little white boxes
  ggplot(aes(x = ICAGR, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Growth rate vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag))
ggsave(paste0("output/risk ICAGR v DD ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

results |>
  ggplot(aes(x = lake, y = bliss, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
            min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Lake ratio vs bliss"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag))
ggsave(paste0("output/risk lake v bliss ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

results |>
  ggplot(aes(x = lake, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Lake ratio vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag))
ggsave(paste0("output/risk lave v DD ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

results |>
  ggplot(aes(x = end_value, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
                   min.segment.length=0, force=0.5, max.iter=10000) +
  scale_x_continuous(labels=scales::dollar_format()) +
  labs(title=paste("Ending value vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag)) 
ggsave(paste0("output/risk end value v DD ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

  # geom_smooth(method = "lm")

results |>
  ggplot(aes(x = ICAGR, y = lake, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.1, "lines"), label.size=0.05, 
       min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Growth rate vs lake ratio"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag)) 
ggsave(paste0("output/risk ICAGR v lake ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)

results |>
  ggplot(aes(x = bliss, y = drawdown, label = lag)) +
  geom_point(shape=4) +
  geom_label_repel(label.padding=unit(0.15, "lines"), label.size=0.05, 
          min.segment.length=0, force=0.5, max.iter=10000) +
  labs(title=paste("Bliss vs drawdowns"),
       subtitle=paste(start_date, "to", end_date, interval, "hr", "EMA:", lag)) 
ggsave(paste0("output/risk bliss v DD ", run_id, run_time, ".pdf"), width=10, height=8, units="in", dpi=300)









