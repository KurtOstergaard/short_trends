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

sum_fun <- function(data, column_name, factor_name) {
  summary_stats <- data |>
    group_by({{ factor_name }}) |>
    summarise(
      mean = mean({{ column_name }}),
      median = median({{ column_name }}),
      sd = sd({{ column_name }}),
      min = min({{ column_name }}),
      max = max({{ column_name }})
    )
  return(summary_stats)
}
#  Good original version
# split_fun <- function(data, column_name, factor_name) {
#   summary_stats <- data |>
#     mutate(category = ifelse({{ column_name }} >= 0, "Positive", "Negative")) |>
#     group_by({{ factor_name }}, category) |>
#     summarise(
#       mean = mean({{ column_name }}),
#       median = median({{ column_name }}),
#       sd = sd({{ column_name }}),
#       min = min({{ column_name }}),
#       max = max({{ column_name }})
#     )
#   return(summary_stats)
# }
  
# ChatGPT new, improved super version
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
# h4s <- read_csv("Aug23-1000.csv", col_names = TRUE)
# a24u <- read_csv("Aug24unix.csv", col_names = TRUE)
# NQ1H <- read_csv("NQ1H.csv", col_names = TRUE)
# NQ1D <- read_csv("NQ1D.csv", col_names = TRUE)
# NQ240 <- read_csv("NQ240.csv", col_names = TRUE)
# h4 <- read_csv("Aug24.csv", col_names = TRUE)

h4 <- read_csv("CME_MINI_NQ1!, 240_48328.csv", col_names = TRUE)
HA_input <- select(h4, time:close)
h4HA <- HAOHLC(HA_input) 
h4HA <- h4HA |> rownames_to_column("time") |>
  mutate(real_open = h4$open,
         real_close = h4$close,
         skid = lead(real_open) - real_close,
         time = as.POSIXlt(time, tz="America/New_York")) |>
  na.omit() |>
  as_tibble() 
  
skid <- 0.5
results <- tibble() 
colnames(results) <- unlist(str_split("j, slow_lag, fast_lag,
            ICAGR, drawdown, bliss, lake, end_val, trade_test", ", "))

# initiate optimization sequence
runs <- expand.grid(lag = seq(2,5,1))
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
    df$drawdown[i] <- df$water[i] / df$equity[i]
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
  }
end_value <- df$equity[nrow(df)] * 20
end_val <- end_value / 1000000
ratio <- end_value/ start_value
start_date <- min(df$time) ; end_date <- max(df$time)
date_range <- as.numeric(difftime(end_date, start_date, units = "days")) / 365.25
ICAGR <- if(ratio <= 0) 0 else log(ratio)/ date_range
drawdown <- max(df$drawdown)
lake <- sum(df$water) / sum(df$equity)
bliss <- ICAGR / drawdown
results[j,1:9] <- as_tibble_row(
  c(j=j, lag=lag, ICAGR=ICAGR, drawdown=drawdown, bliss=bliss, lake=lake, 
    end_value=end_value, end_val=end_val, trade_test=trade_test),
  .name_repair = "universal")

}       # optimization loop end
end_time <- Sys.time() ;forever <- end_time - start_time
secs <- forever  / nrow(runs)
sprintf("Yo, %1.2f total time and %1.4f per run, %i runs", forever, secs, nrow(runs))




zz <- split_fun(trades, trade_pnl)
zz

# the promised land of pretty graphs
df <- df |>
  mutate(date = as.POSIXct(time)) 
df |>
  ggplot(aes(x = date, y = real_close)) +
  geom_line(size = 1, alpha = 0.5)

ggplot() +
  geom_point(data=df, aes(x=date, y=close), alpha=0.2) +
  geom_segment(data=trades, aes(x=buy_date, y=buy_price, xend=sell_date,
                                yend=sell_price, size = 1, color="black"))

df |>        # lake
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin=equity, ymax=highwater, x=date, fill = "band"), alpha = 0.9)+
  scale_color_manual("", values="grey12")+
  scale_fill_manual("", values="red")

df |>       # equity and highwater
  mutate(close_big = close*10500) |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = equity),size = 1,  alpha = 0.9) +
  geom_point(aes(y=highwater), size=1, shape = 4, alpha=0.1, color="black") +
  geom_line(aes(y=close_big), size=1, alpha=0.2)

# check for second scale
# https://stackoverflow.com/questions/59956953/plotting-secondary-axis-using-ggplot

df |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = equity), size = 1, alpha = 0.8) +
  geom_line(aes(y = highwater), size=1, alpha=0.2)

rzlt <- results
 

rzlt |>
  ggplot(aes(x = ICAGR, y = drawdown)) +
  geom_point(size = 3, shape = 4)

rzlt |>
  ggplot(aes(x = lake, y = bliss)) +
  geom_point(size = 3, shape = 4)

rzlt |>
  ggplot(aes(x = lake, y = drawdown)) +
  geom_point(size = 3, shape = 4)

#  geom_smooth(method = "lm")

rzlt |>
  ggplot(aes(x = ICAGR, y = lake)) +
  geom_point(size = 3, shape = 4)

rzlt |>
  ggplot(aes(x = ICAGR, y = bliss)) +
  geom_point(size = 3, shape = 4)


path <- getwd()
now <- Sys.time()
filename <- paste(path, "/results ", now, ".csv", sep="")
write_csv(results, filename)

################ New Graph idea

# Sample market price data (replace with your actual market data)
# set.seed(42)
# start_date <- as.POSIXct("2023-01-01", tz = "UTC")
# time_intervals <- seq(start_date, length.out = 10 * 365 * 6, by = "4 hours")
# market_prices <- data.frame(
#   time = time_intervals,
#   price = cumsum(runif(length(time_intervals), min = -0.5, max = 0.5))
# )

# Sample trade data (replace with your actual trade data)
# trades <- data.frame(
#   buy_date = sample(time_intervals, size = 10),
#   buy_price = runif(10, min(market_prices$price), max(market_prices$price)),
#   sell_date = sample(time_intervals, size = 10),
#   sell_price = runif(10, min(market_prices$price), max(market_prices$price))
# )

# Define a function to create the ggplot graph
create_trade_plot <- function(df, trades) {
  trades <- trades |>
    mutate(
      buy_date_day_before = buy_date - days(1),
      sell_date_day_after = sell_date + days(1),
      trade_index = row_number()  # Create an index variable for trade facets
    )
  
  p <- ggplot() +
    geom_line(data = df, aes(x = date, y = real_close)) +
    geom_segment(data = trades, aes(x = buy_date, xend = sell_date, y = buy_price, yend = sell_price),
                 color = "red", size = 1) +
    facet_wrap(~ trade_index, scales = "free_x", ncol = 1) +
    labs(x = "Time", y = "Price") +
    theme_minimal()
  
  return(p)
}

# Create the ggplot graph
plot <- create_trade_plot(df, trades)

# Print the graph
print(plot)

######################## End of new graph idea






















###### stop here
pnl <- rzlt |>
  select(slow_lag, fast_lag, end_val) |>
  pivot_wider(names_from = slow_lag, values_from = end_val) |>
  as.matrix()
blip <-  sort(unique(pnl[,1]))
rownames(pnl) <- blip
pnl <- pnl[,-1]


plot_ly(z = ~pnl) |> add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)

happy <- rzlt |>
  select(slow_lag, fast_lag, bliss) |>
  pivot_wider(names_from = slow_lag, values_from = bliss) |>
  as.matrix()
happy <- happy[,-1]
xax <- sort(unique(rzlt$slow_lag))
yax <- sort(unique(rzlt$fast_lag))

plot_ly(z = happy, x=xax, y=yax, type = "surface") |>
  layout(scene = list(
    xaxis=list(title="slow"),
    yaxis=list(title="fast", autorange="reversed"),
    zaxis=list(title="Bliss")
  ))


fig <- plot_ly(z = happy, x=xax, y=yax) |> add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)
fig <- fig |>  layout(
  scene = list(
    xaxis=list(title="slow"),
    yaxis=list(title="fast", autorange="reversed"),
    zaxis=list(title="Bliss")
  ))

fig


plot_ly(z = ~happy) |> add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)

cc <-  cor(rzlt$ICAGR, rzlt$drawdown, method="pearson")
cc
plot_ly(z = volcano, type = "surface")

# https://win-vector.com/2015/07/27/efficient-accumulation-in-r/
# notes on efficient accumulation

trades |>       # trades
  ggplot(aes(x = date)) +
  geom_segment(aes(x=buy_date, y=buy_price, xend=sell_date, yend=sell_price,
                   color="black"))

trades |>       # trades
  ggplot(aes(x = date)) +
  geom_segment(aes(x=buy_date, y=buy_price, xend=sell_date, yend=sell_price,
                   color="black"))




















sketch <- tibble() |>
  colnames(unlist(list("m30", "h1", "h2", "h3", "h4", "h6")))

evry <- bind_rows(h6, h4, h3, h2, h1, m30, .id="id") 

evry1 <- evry |>
  select(time:close) |>
  distinct(.keep_all = FALSE) |>
  arrange(time)


test <- h6 |>
  filter(time >"2023-08-15" & time < "2023-08-17")









