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

#################### Heikin Ashi
# https://stackoverflow.com/questions/53941102/quantmod-heikin-ashi-plotting-unavailable
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

# p <- getSymbols("SPY")
# microbenchmark(HAOHLC(p), times = 100L, unit = "ms")
# microbenchmark(HAOHLC(p), heikin_ashi(p), times = 100L, unit = "ms")

# Unit: milliseconds
# expr       min          lq        mean     median         uq        max neval
# HAOHLC(p)   0.36409   0.4275205   0.5198086   0.502614   0.552392   1.378134   100
# heikin_ashi(p) 563.33925 582.6144955 609.0082902 591.338550 620.179235 802.885348   100

################## Heikin Ashi 2 


heikin_ashi <- function(data) {
  
  if(!quantmod::is.OHLC(data)) stop("data must contain OHLC columns")
  
  heikin_close <- xts::xts(Matrix::rowMeans(quantmod::OHLC(data)), order.by = index(data))
  heikin_open  <- quantmod::Op(data)
  
  # need a loop: heiki ashi open is dependent on the previous value
  for(i in 2:nrow(data)) {
    heikin_open[i] <- (heikin_open[i-1] + heikin_close[i-1]) / 2
  }
  
  heikin_high <- xts::xts(apply(cbind(quantmod::Hi(data), heikin_open, heikin_close), 1, max), order.by = index(data))
  heikin_low <- xts::xts(apply(cbind(quantmod::Lo(data), heikin_open, heikin_close), 1, min), order.by = index(data))
  
  out <- merge(heikin_open, heikin_high, heikin_low, heikin_close)
  out <- setNames(out, c("Open", "High", "Low", "Close"))
}



ADM <- getSymbols("ADM", from = "2018-10-01", auto.assign = FALSE)
ha_ADM <- heikin_ashi(AMD)
chartSeries(ha_ADM) # or chart_Series or rtsplot or geom_candlestick

############################# end of Heikin Ashi



m30 <- read_csv("CME_MINI_NQ1!, 30_ac16d.csv", col_names = TRUE)
h1 <- read_csv("CME_MINI_NQ1!, 60_f371d.csv", col_names = TRUE)
h2 <- read_csv("CME_MINI_NQ1!, 120_925b9.csv", col_names = TRUE)
h3 <- read_csv("CME_MINI_NQ1!, 180_62ca4.csv", col_names = TRUE)
h4 <- read_csv("CME_MINI_NQ1!, 240_48328.csv", col_names = TRUE)
h6 <- read_csv("CME_MINI_NQ1!, 360_2a174.csv", col_names = TRUE)
spec(h6)

h4HA <- HAOHLC(h4)
h4HA <- h4HA |> rownames_to_column("time") 
h4HA$time <- as.POSIXct(h4HA$time, tz=Sys.timezone())

h4b <- h4 
h4c <- full_join(h4b, h4HA)
  
  
  
sketch <- tibble() |>
  colnames(unlist(list("m30", "h1", "h2", "h3", "h4", "h6")))

test1 <- filter(h1, time > "2023-08-13" & time < "2023-08-18")
test6 <- filter(h6, time > "2023-08-13" & time < "2023-08-18")



evry <- bind_rows(h6, h4, h3, h2, h1, m30, .id="id") 

evry1 <- evry |>
  select(time:close) |>
  distinct(.keep_all = FALSE) |>
  arrange(time)


test <- h6 |>
  filter(time >"2023-08-15" & time < "2023-08-17")









