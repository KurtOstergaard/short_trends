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
# library(clock)

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

# Start dates: 1hr starts 1/1/20, 2hr 1/2/17, 3hr & 4hr 8/14/13, 6hr 1/2/04
# m30 <- read_csv("CME_MINI_NQ1!, 30_ac16d.csv", col_names = TRUE)
# h1 <- read_csv("CME_MINI_NQ1!, 60_f371d.csv", col_names = TRUE)
# h2 <- read_csv("CME_MINI_NQ1!, 120_925b9.csv", col_names = TRUE)
# h3 <- read_csv("CME_MINI_NQ1!, 180_62ca4.csv", col_names = TRUE)
# h4 <- read_csv("CME_MINI_NQ1!, 240_48328.csv", col_names = TRUE)
# h6 <- read_csv("CME_MINI_NQ1!, 360_2a174.csv", col_names = TRUE)
# spec(h6)
h4 <- read_csv("Aug23-1000.csv", col_names = TRUE)

h4HA <- HAOHLC(h4) 
h4HA <- h4HA |> rownames_to_column("time") |>
  mutate(real_open = h4$open,
         real_close = h4$close,
         skid = lead(real_open) - real_close,
         time = as.POSIXlt(time, tz="America/New_York")) |>
  na.omit() |>
  as_tibble() 
  
# h4NYC <- with_tz(h4, tz = "America/New_York")
# time = as.POSIXct(time, tz=Sys.timezone())

# max(h4HA$spread) # 171 Wow!
# min(h4HA$spread) # -287.5 double Wow!
# holy fuck! Make a histogram of this craziness to estimate buy skid
skid <- 1

results <- tibble() 
colnames(results) <- unlist(str_split("j, slow_lag, fast_lag,
            ICAGR, drawdown, bliss, lake, end_val, trade_test", ", "))

# Scenario or Optimization?
j <- 18     # Scenario run
# for (j in seq(3, 5, 1)) {   #initiate optimization sequence
  df <- h4HA
  df$lag <- ewmaRcpp(df$real_close, j)
  
  df<- df |>                # create trade signal 
    mutate(cross = ha.close - lag,
           on = ifelse(cross > 0 & lag(cross) < 0, 1, 0),
           off = ifelse(cross < 0 & lag(cross) > 0, -1, 0))
  df <- slice(df, 30:n()) # drop first 30 rows to warm up EMA & Heikin Ashi
  
  if(df$cross[1] > 0) df$on[1] <- 1  # first row
  df <- df|>                         # trade signal details 
    mutate(signal = on + off,
          buy_date = ifelse(on == 1, as_datetime(lead(time)), 0),
          buy_price = ifelse(on == 1, real_close + skid, 0),
          buy_amount = 0,
          sell_date = ifelse(off == -1, as_datetime(lead(time)), 0),
          sell_price = ifelse(off == -1, real_close - skid, 0))

  if (df$cross[nrow(df)] > 0) {    # close out trade if long at EOF
  df$off[nrow(df)] <- -1
  df$signal[nrow(df)] <- -1
  df$sell_date[nrow(df)] <- df$time[nrow(df)]
  df$sell_price[nrow(df)] <- df$real_close[nrow(df)] - skid
}
 
  start_value <- 1e5 ; df$closed_pnl <- 1e5 ; df$highwater <- 1e5
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
  # summary trade table extract (or extrusion?)
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
  }
  z <- df$equity[nrow(df)]
z
zz <- sum(trades$trade_pnl)
zz
# calc win/lose count and avg win/lose $        Start here!
ntrades <- nrow(trades)
count_wins <- trades |>
  filter(trades$trade_pnl > 0) |>
  nrow()
count_lose <- trades |>
  filter(trades$trade_pnl <= 0) |>
  nrow()

win_pnl <- trades |>
  filter(trades$trade_pnl > 0) |>
  sum()


trades |>
  select(trades$trade_pnl > 0) |>
  win_count <- nrow |>
  win_pnl = sum()

# }   # optimization end














sketch <- tibble() |>
  colnames(unlist(list("m30", "h1", "h2", "h3", "h4", "h6")))

evry <- bind_rows(h6, h4, h3, h2, h1, m30, .id="id") 

evry1 <- evry |>
  select(time:close) |>
  distinct(.keep_all = FALSE) |>
  arrange(time)


test <- h6 |>
  filter(time >"2023-08-15" & time < "2023-08-17")









