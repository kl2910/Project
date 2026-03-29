install.packages("tidyquant")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("GGally")
install.packages("quadprog")

#
library(tidyquant)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GGally)
library(quadprog)

# Define tickers (your 10 stocks)
tickers <- c("COST","AMZN","AU","FET","AAPL","WMT","GOOG","NVDA","PEP","MU")

# IMPORTANT: pick a start date that works for all tickers (avoid too many NAs)
start_date <- "2018-01-01"

# Get monthly returns for all stocks
returns_data <- tq_get(tickers,
                       get = "stock.prices",
                       from = start_date) %>%
  group_by(symbol) %>%
  tq_transmute(
    select     = adjusted,
    mutate_fun = periodReturn,
    period     = "monthly",
    col_rename = "monthly_return"
  ) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = symbol, values_from = monthly_return)

# Remove rows with missing data
returns_data_clean <- na.omit(returns_data)

# (optional) print correlation matrix
cor_mat <- cor(returns_data_clean[, -1])
print(round(cor_mat, 3))

# View a correlation scatterplot matrix
ggpairs(returns_data_clean[, -1],
        title = "Monthly Return Correlation Matrix: 10 Stocks")

# Market proxy: S&P 500 (^GSPC)
sp500 <- tq_get("^GSPC", get = "stock.prices", from = start_date) %>%
  tq_transmute(
    select     = adjusted,
    mutate_fun = periodReturn,
    period     = "monthly",
    col_rename = "market_return"
  )

# Risk-free rate (monthly)
rf_annual  <- 0.0432
rf_monthly <- rf_annual / 12

# Function to compute Sharpe + Jensen alpha for ONE ticker
calc_metrics <- function(ticker) {
  
  stock <- tq_get(ticker, get = "stock.prices", from = start_date) %>%
    tq_transmute(
      select     = adjusted,
      mutate_fun = periodReturn,
      period     = "monthly",
      col_rename = "stock_return"
    )
  
  returns <- left_join(stock, sp500, by = "date") %>%
    na.omit() %>%
    mutate(
      excess_stock  = stock_return  - rf_monthly,
      excess_market = market_return - rf_monthly
    )
  
  sharpe <- mean(returns$excess_stock, na.rm = TRUE) /
    sd(returns$stock_return,  na.rm = TRUE)
  
  model <- lm(excess_stock ~ excess_market, data = returns)
  alpha <- coef(model)[1]
  beta  <- coef(model)[2]
  
  tibble(
    ticker = ticker,
    sharpe_monthly = sharpe,
    alpha_monthly  = alpha,
    beta           = beta,
    alpha_annualized = (1 + alpha)^12 - 1
  )
}

# Compute for all 10 stocks
metrics_table <- map_dfr(tickers, calc_metrics)

# Print nicely
metrics_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()

# Use the cleaned wide returns matrix (exclude date)
R <- as.matrix(returns_data_clean %>% select(-date))

mu    <- colMeans(R)
Sigma <- cov(R)
rf    <- rf_monthly

# Negative Sharpe (optim minimizes)
neg_sharpe <- function(w, mu, Sigma, rf) {
  w <- w / sum(w)
  port_ret <- sum(w * mu)
  port_sd  <- sqrt(t(w) %*% Sigma %*% w)
  - (port_ret - rf) / port_sd
}

n_assets <- ncol(R)
w0 <- rep(1 / n_assets, n_assets)

opt <- optim(
  par    = w0,
  fn     = neg_sharpe,
  mu     = mu,
  Sigma  = Sigma,
  rf     = rf,
  method = "L-BFGS-B",
  lower  = rep(0, n_assets),
  upper  = rep(1, n_assets)
)

weights <- opt$par / sum(opt$par)
names(weights) <- colnames(R)

cat("Optimal Portfolio Weights (Max Sharpe Ratio):\n")
print(round(weights, 4))

portfolio_return <- sum(weights * mu)
portfolio_sd     <- sqrt(t(weights) %*% Sigma %*% weights)
portfolio_sharpe <- (portfolio_return - rf) / portfolio_sd

cat("\nPortfolio Expected Monthly Return:", round(portfolio_return, 4), "\n")
cat("Portfolio Monthly Volatility:", round(portfolio_sd, 4), "\n")
cat("Portfolio Sharpe Ratio:", round(portfolio_sharpe, 4), "\n")


       