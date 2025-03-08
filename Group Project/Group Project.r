# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Asset Pricing and Portfolio Management | GROUP PROJECT  | Afonso Casaonva, 20240795 | Francisco Perestrello, 20241560 | João Pardal, 20240796 | Nuno Vieira, 20241111
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls(all.names = TRUE))
graphics.off()
close.screen(all.screens = TRUE)
erase.screen()
windows.options(record=TRUE)

library(pacman)
if (!require("pacman")) {install.packages("pacman")}

p_load(xts)                   # to manipulate time series of stock data
p_load(quantmod)              # to download stock data
p_load(PerformanceAnalytics)  # to compute performance measures
p_load(tseries)               # time series analysis
p_load(FinTS)                 # to compute an ARCH test
p_load(dplyr)                 # for data manipulation
p_load(rcompanion)            # for statistical analysis
p_load(GLDEX)                 # for moments (mean, variance, skewness, kurtosis)
p_load(rugarch)               # to use GARCH Models
p_load(CVXR)                  # Convex Optimization in R
p_load(riskParityPortfolio)   # RPP
p_load(ggplot2)               # plot
p_load(portfolioBacktest)     # portfolio Backtesting

# -------------------------------------------------------------------------------------------
#                                       2. Group Part
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# Load stock price data and compute returns
# -------------------------------------------------------------------------------------------

# Download data from Yahoo Finance
start_date <- "2010-01-01"
end_date   <- "2024-10-31"

Tickers <- c("SPY","GOOGL","NFLX","AAPL","AMZN","META","GLD","NVDA","IBM","TXN","ASML","DECK","V","MA","SNPS","JPM","AVGO","AMAT","TLT","NVO","LLY","AXON","XYL","GBTC","NTES","NEE","JNJ","PG","LMT","TSLA","LRCX")
# S&P 500, Bitcoin, Gold, 20+Yr US Bond and 27 Stocks ranging from Healthcare, to Financials, to Consumer electronic, Semiconductors, Software and Industrials

prices_d <- xts()
for (i in 1:length(Tickers)) {
  tmp <- Ad(getSymbols(Tickers[i], from = start_date, to = end_date, 
                       periodicity = "daily", auto.assign = FALSE))  # Daily adjusted prices
  tmp <- na.approx(tmp, na.rm = FALSE)  # Interpolate NAs
  prices_d <- cbind(prices_d, tmp)
  prices_d <- na.omit(prices_d) # Remove NAs that last even after interpolation
}

colnames(prices_d) <- Tickers
tclass(prices_d) <- "Date"

# Calculate weekly and monthly prices
prices_w <- to.period(prices_d, period = "weeks", k = 1, OHLC = FALSE) # The method uses Last Trading Day of Week
prices_m <- to.period(prices_d, period = "months", k = 1, OHLC = FALSE) # The method uses Last Trading Day of Month

# Compute returns
X_log_d <- CalculateReturns(prices_d, "log")[-1]
X_log_w <- CalculateReturns(prices_w, "log")[-1]
X_log_m <- CalculateReturns(prices_m, "log")[-1]

# Annualized returns, standard deviation
t(table.AnnualizedReturns(X_log_d, scale = 252, Rf = 0.0, geometric = T))
t(table.AnnualizedReturns(X_log_w, scale = 52, Rf = 0.0, geometric = T))
t(table.AnnualizedReturns(X_log_m, scale = 12, Rf = 0.0, geometric = T))

# -------------------------------------------------------------------------------------------
# Investigate for the absence of auto-correlation
# -------------------------------------------------------------------------------------------

# Perform the Ljung-Box autocorrelation test
ljung_box_test <- function(data, lag=30) {
  lb_test <- Box.test(data, lag = lag, type = "Ljung-Box")
  return(lb_test$p.value)
}

autocorrelation <- data.frame()
for (ticker in Tickers) {
  p_value_d <- ljung_box_test(X_log_d[, ticker])
  p_value_w <- ljung_box_test(X_log_w[, ticker])
  p_value_m <- ljung_box_test(X_log_m[, ticker])
  autocorrelation <- rbind(autocorrelation, 
                           data.frame(daily = round(p_value_d, 4), 
                                      weekly = round(p_value_w, 4), 
                                      monthly = round(p_value_m, 4)))
}
row.names(autocorrelation) <- Tickers

autocorrelation # null hypothesis being that the data exhibits no autocorrelation at any lag up to the specified lag. -> p-value below 0.05 rejects the null


############ Relevant plotting code ############ 

# Plotting autocorrelation functions
plot_autocorrelation <- function(data, title) {
  par(mfrow = c(2,1))
  acf(data, main = paste("ACF -", title))
  pacf(data, main = paste("PACF -", title))
}

# Example usage
plot_autocorrelation(X_log_d$SPY, "SPY Log Returns")

# -------------------------------------------------------------------------------------------
# Investigate for Fat Tails
# -------------------------------------------------------------------------------------------

jarque_bera_test <- function(data) {
  jb_test <- jarque.bera.test(data)
  return(jb_test$p.value)
}

non_normality <- data.frame()
for (ticker in Tickers) {
  p_value_d <- jarque_bera_test(X_log_d[, ticker])
  p_value_w <- jarque_bera_test(X_log_w[, ticker])
  p_value_m <- jarque_bera_test(X_log_m[, ticker])
  non_normality <- rbind(non_normality, 
                         data.frame(daily = round(p_value_d, 4), 
                                    weekly = round(p_value_w, 4), 
                                    monthly = round(p_value_m, 4)))
}
row.names(non_normality) <- Tickers

non_normality # null hypothesis being that the data exhibits normality -> p-value below 0.05 rejects the null


############ Relevant plotting code ############

# Graphical comparison of the empirical quantiles of a data sample to those from a specified reference distribution
par(mfrow=c(2,1))
chart.QQPlot(X_log_w$GLD, main = "GLD",
             line=c("quartiles"), distribution = 'norm', envelope=0.95)
chart.QQPlot(X_log_m$GLD, main = "GLD",
             line=c("quartiles"), distribution = 'norm', envelope=0.95)


# Normal Distribution overlay on Histogram
par(mfrow=c(1,1))
plotNormalHistogram(X_log_d$GBTC, prob = TRUE,
                    breaks=50, col='cornflowerblue',
                    linecol='red', main = "GE",
                    length = 10000, xlab='Returns')

# -------------------------------------------------------------------------------------------
# Investigate for Skewness
# -------------------------------------------------------------------------------------------

mom_d = apply(X_log_d, 2, fun.moments.r, normalise = "N")
rbind(mom_d, Excess.Kurtosis=mom_d[4,]-3)

mom_w = apply(X_log_w, 2, fun.moments.r, normalise = "N") 
rbind(mom_w, Excess.Kurtosis=mom_w[4,]-3)

mom_m = apply(X_log_m, 2, fun.moments.r, normalise = "N")
rbind(mom_m, Excess.Kurtosis=mom_m[4,]-3)


############ Relevant plotting code ############

# Function to plot returns density function
plot_returns_density <- function(returns, title = "Returns Density", col = "cornflowerblue") {
  par(mfrow=c(1,1))
  # Create density plot
  plot(density(returns, na.rm = TRUE), main = title, xlab = "Returns", col = col, lwd = 2)
  # Add a vertical line at x=mean for reference
  abline(v = mean(returns), col = "red", lty = 2)
}

# Example usage:
plot_returns_density(X_log_w$GOOGL)


# Function to create boxplot of returns
plot_returns_boxplot <- function(returns, title = "Returns Distribution") {
  par(mfrow=c(1,1))
  # Create boxplot
  boxplot(returns, 
          main = title,
          ylab = "Returns",
          col = "cornflowerblue",
          border = "darkblue",
          notch = TRUE,
          outcol = "red",
          outpch = 20)
  
  # Add horizontal line at y=0 for reference
  abline(h = 0, col = "red", lty = 2)
}

# Example usage:
plot_returns_boxplot(X_log_w$SPY, "Google Daily Returns Distribution")


# Function to plot histogram
plot_returns_histogram <- function(returns, title = "Returns Histogram", breaks = 50, col = "cornflowerblue") {
  hist(returns, main = title, col = col, xlab = "Returns", breaks = breaks)
  # Add horizontal line at y=0 for reference
  abline(v = 0, col = "red", lty = 2)
}

# Example usage:
plot_returns_histogram(X_log_d$SPY)


# -------------------------------------------------------------------------------------------
# Investigate for Volatility Clustering
# -------------------------------------------------------------------------------------------

# The Autoregressive Conditional Heteroskedasticity (ARCH) test can test for volatility clustering by checking if there’s conditional heteroscedasticity (time-varying volatility)
arch_test <- function(data) {
  arch_test <- ArchTest(data, lags = 9)
  return(arch_test$p.value)
}

volatility_clustering <- data.frame()
for (ticker in Tickers) {
  p_value_d <- arch_test(X_log_d[, ticker])
  p_value_w <- arch_test(X_log_w[, ticker])
  p_value_m <- arch_test(X_log_m[, ticker])
  volatility_clustering <- rbind(volatility_clustering, 
                                data.frame(daily = round(p_value_d, 4), 
                                           weekly = round(p_value_w, 4), 
                                           monthly = round(p_value_m, 4)))
}
row.names(volatility_clustering) <- Tickers

volatility_clustering # null hypothesis being that the data exhibits no volatility clustering -> p-value below 0.05 rejects the null


############ Relevant plotting code ############

plot_volatility_cluster <- function(data, title) {
  rolling_volatility <- rollapply(data, width = 14, FUN = sd, align="right", fill = NA)
  par(mfrow=c(1,1))
  plot(data, main=title, lwd=1, col="blue") # Plot the returns
  lines(rolling_volatility, main = "20-Day Rolling Volatility", lwd=1, col = "red") # Plot the rolling volatility
}

plot_volatility_cluster(X_log_m$GOOGL, "Google Log Returns") # Example

# -------------------------------------------------------------------------------------------
# Investigate for Leverage Effects
# -------------------------------------------------------------------------------------------

leverage_effect <- function(data) {
  df <- data.frame(log_returns = data)
  df <- df %>% mutate(lagged_returns = lag(data), squared_returns = data^2)
  correlation <- cor(df$lagged_returns, df$squared_returns, use = 'complete.obs') # use=complete.obs handles NAs by casewise deletion
  return(correlation[1,1])
}

leverage <- data.frame()
for (ticker in Tickers) {
  corr_d <- leverage_effect(X_log_d[, ticker])
  corr_w <- leverage_effect(X_log_w[, ticker])
  corr_m <- leverage_effect(X_log_m[, ticker])
  leverage <- rbind(leverage,
                    data.frame(daily = round(corr_d, 4),
                               weekly = round(corr_w, 4),
                               monthly = round(corr_m, 4)))
}
row.names(leverage) <- Tickers

leverage # a negative correlation suggests that negative returns are associated with higher future volatility, indicating the presence of a leverage effect


############ Relevant plotting code ############

news_impact_curve <- function(data, ticker) {
  # Define an EGARCH model specification
  spec <- ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm")

  # Fit the EGARCH model to log-returns
  egarch_model <- ugarchfit(spec = spec, data = data[,ticker])

  # Extract model parameters
  alpha <- coef(egarch_model)["alpha1"]  # asymmetric term
  beta <- coef(egarch_model)["beta1"]
  omega <- coef(egarch_model)["omega"]
  gamma <- coef(egarch_model)["gamma1"]  # leverage effect term

  # Create a range of return values to simulate the news impact
  news <- seq(-0.05, 0.05, length.out = 100)

  # Calculate conditional volatility response for each news level
  NIC <- exp(omega + (alpha * abs(news)) + (gamma * news) + beta)

  # Plot the News Impact Curve
  par(mfrow=c(1,1))
  plot(news, NIC, type = "l", col = "blue", lwd = 2,
      main = "News Impact Curve",
      xlab = "News (Return Shock)",
      ylab = "Conditional Volatility Response")
  abline(v = 0, col = "red", lty = 2)
}

# Example usage
news_impact_curve(X_log_m, 'JNJ')

# -------------------------------------------------------------------------------------------
# Investigate for Conditional Non-Normality
# -------------------------------------------------------------------------------------------

# Split data into yearly periods and test for normality
conditional_non_normality <- function(log_returns) {
  normality_results <- data.frame()
  for (ticker in Tickers) {
    # Extract returns for current ticker
    returns <- log_returns[, ticker]
    
    # Get years in the data
    years <- unique(format(index(returns), "%Y"))
    
    # Test each year
    for (year in years) {
      # Get data for current year
      yearly_data <- returns[format(index(returns), "%Y") == year]
      
      # Perform Shapiro-Wilk test
      test_result <- shapiro.test(coredata(yearly_data)) # coredata retrieves the numeric vector of the xts object
      
      # Store results
      normality_results <- rbind(normality_results,
                                data.frame(Ticker = ticker,
                                        Year = year,
                                        P_Value = round(test_result$p.value, 4),
                                        Is_Normal = test_result$p.value > 0.05))
    }
  }
  return(normality_results)
}

normality_results_d <- conditional_non_normality(X_log_d) # null hypothesis being that the data exhibits normality -> p-value below 0.05 rejects the null
normality_results_w <- conditional_non_normality(X_log_w)
normality_results_m <- conditional_non_normality(X_log_m)


# Calculate percentage of normal vs non-normal distributions for the daily example
total_tests <- nrow(normality_results_d)
normal_count <- sum(normality_results_d$Is_Normal)
non_normal_count <- total_tests - normal_count

normal_count; round(100 * normal_count/total_tests, 2)
non_normal_count; round(100 * non_normal_count/total_tests, 2)


############ Relevant plotting code ############

plot_non_normality <- function(data, normality_results, ticker, freq) {
  # Convert xts to data frame for ggplot
  df <- data.frame(
    Date = index(data[,ticker]),
    Returns = as.numeric(data[,ticker])
  )
  # Filter normality results for this ticker
  results <- normality_results[normality_results$Ticker == ticker,]
  # Create year ranges for shading
  year_ranges <- data.frame(
    start = as.Date(paste0(results$Year, "-01-01")),
    end = as.Date(paste0(results$Year, "-12-31")),
    Is_Normal = results$Is_Normal
  )
  # Create plot
  ggplot(df, aes(x = Date, y = Returns)) +
    geom_rect(data = subset(year_ranges, !Is_Normal),
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf,
              fill = "red", alpha = 0.2,
              inherit.aes = FALSE) +
    geom_line(color = "blue", linewidth = 0.5) +
    labs(title = paste("Returns Distribution for", ticker, "(", freq, ")"),
         subtitle = "Red shading indicates non-normal periods",
         x = "Date", 
         y = "Log Returns")
}

# Example usage
plot_daily <- plot_non_normality(X_log_d, normality_results_d, "GBTC", "daily") 
plot_weekly <- plot_non_normality(X_log_w, normality_results_w, "GBTC", "weekly")
plot_monthly <- plot_non_normality(X_log_m, normality_results_m, "GBTC", "monthly")

p_load(gridExtra)
grid.arrange(plot_daily, plot_weekly, plot_monthly, nrow = 1, ncol = 3)

# -------------------------------------------------------------------------------------------
#                                       2. Backtesting
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# Portfolio definitions
# -------------------------------------------------------------------------------------------

# Note from portfolioBacktest() documentation:
# Each portfolio design is easily defined as a function that takes as input a window of the stock prices and outputs the portfolio weights.

# Equally Weighted Portfolio (EWP)
EWP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  N <- ncol(X)
  w <- rep(1/N, N)
  names(w) <- colnames(X)
  return(w)
}

# Markowitz Mean-Variance Portfolio (MVP)
MVP_portfolio_fun <- function(dataset, lambda=0.5, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu    <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)       # compute the SCM
  # design mean-variance portfolio
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lambda*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- colnames(Sigma)
  return(w)
}

# Global Minimum Variance Portfolio (GMVP)
GMVP_portfolio_fun <- function(dataset, Sigma = NULL, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  # If Sigma is not provided, compute it from log_returns                  # Additional flexibility to allow for MDCP usage
  if (is.null(Sigma)) {
    Sigma <- cov(X)
  }
  w <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- colnames(Sigma)
  return(w)
}

# Maximum Sharpe Ratio Portfolio (MSRP) (in convex form)
MSRP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu <- colMeans(X)
  Sigma <- cov(X)
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1))     # Minimize varP(w) s.t. E(Rp)=1, w=w/sum(w)
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w_)/sum(result$getValue(w_)))
  names(w) <- colnames(Sigma)
  return(w)
}

# Inverse Volatility Portfolio (IVP)
IVP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)
  sigma <- sqrt(diag(Sigma))
  w <- 1/sigma
  w <- w/sum(w)
  return(w)
}

# Vanilla Risk Parity Portfolio (RPP)
RPP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)
  rpp <- riskParityPortfolio(Sigma)  # Returns weights and relative risk contributions
  return(rpp$w)
}

# Most Diversified Portfolio (MDP)
MDP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)
  vol <- sqrt(diag(Sigma))
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(vol) %*% w_ == 1))     # Minimize varP(w) s.t. E(σ)=1, w=w/sum(w)
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w_)/sum(result$getValue(w_)))
  names(w) <- colnames(Sigma)
  return(w)
}

# Maximum Decorrelation Portfolio (MDCP)
MDCP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma = cov(X)
  C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
  colnames(C) <- colnames(Sigma)
  return(GMVP_portfolio_fun(dataset, Sigma = C))
}


# -------------------------------------------------------------------------------------------
# Backtesting
# -------------------------------------------------------------------------------------------

# Download prices again because financialDataResample() needs that the data argument matches the structure returned by the function stockDataDownload() =/= getSymbols()
prices <- stockDataDownload(stock_symbols = Tickers,
                            index_symbol = "^GSPC",
                            from = start_date, to = end_date)


save(prices, file = "stockdata_from_2010-01-01_to_2024-10-31.Rdata")
#load('stockdata_from_2010-01-01_to_2024-10-31.Rdata')    # Once loaded the data for the first time, save time by running this line of code instead

# Resample 100 times, each with 10 stocks and 2-year consecutive data
set.seed(123)

my_dataset_list <- financialDataResample(prices, 
                                         N_sample = 10,     # Desired number of financial instruments in each resample
                                         T_sample = 252*2,  # Desired length of each resample 
                                         num_datasets = 100) # Number of resampled datasets

portfolios <- list("EWP" = EWP_portfolio_fun,
                   "GMVP" = GMVP_portfolio_fun,
                   "MVP" = MVP_portfolio_fun,
                   "MSRP" = MSRP_portfolio_fun,
                   "IVP" = IVP_portfolio_fun,
                   "RPP" = RPP_portfolio_fun,
                   "MDP" = MDP_portfolio_fun,
                   "MDCP" = MDCP_portfolio_fun)

bt <- portfolioBacktest(portfolios, my_dataset_list,
                        benchmark = ("index"),
                        rebalance_every = 252*0.5,           # arbitrary
                        optimize_every = 252*0.5,            # arbitrary
                        lookback = 252*0.5,                  # arbitrary
                        shortselling = F,
                        cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4)) # arbitrary

# Summary of portfolios performance measures
res_sum <- backtestSummary(bt, summary_fun = median, show_benchmark = TRUE) # Show median values on the table
res_sum <- backtestSummary(bt, summary_fun = function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE), show_benchmark = TRUE) # Show max values on the table

summaryTable(res_sum, type = "DT", 
             order_col = "Sharpe ratio",   # performance measure to be used to sort the rows
             order_dir = "desc")

# BoxPlot of portfolio simulations' Sharpe Ratios:
backtestBoxPlot(bt[c("MVP","GMVP","MSRP","index")], measure = "Sharpe ratio") # This example selects specific portfolios

# BoxPlot of portfolio simulations' Max drawdowns:
backtestBoxPlot(bt, measure = "max drawdown")                                 # This example selects all portfolios

# Cumulative Return or wealth plot of a single backtest
backtestChartCumReturn(bt, dataset_num = 1)

# Drawdown plot of a single backtest
backtestChartDrawdown(bt, dataset_num = 1)

# Portfolio allocation evolution of a particular portfolio over a single backtest
backtestChartStackedBar(bt, "MDCP", dataset_num = 1, legend = TRUE)


# Function to plot different portfolios on the efficient frontier
plot_efficient_frontier_many <- function(dataset, weights_list, portfolios) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu <- colMeans(log_returns)  # compute mean vector
  Sigma <- cov(log_returns)  # compute covariance matrix
  
  # Create a range of risk aversion parameters
  lambdas <- seq(0, 100, length.out = 200)
  
  # Initialize vectors to store portfolio metrics
  portfolio_returns <- numeric(length(lambdas))
  portfolio_risks <- numeric(length(lambdas))
  
  # Calculate portfolios for different risk aversion parameters
  for(i in seq_along(lambdas)) {
    # Get the optimized weights
    w_opt <- MVP_portfolio_fun(dataset, lambda = lambdas[i]) # Use the mean-variance portfolio to draw the frontier
    
    # Compute portfolio return and risk
    portfolio_returns[i] <- t(mu) %*% w_opt
    portfolio_risks[i] <- sqrt(t(w_opt) %*% Sigma %*% w_opt)
  }
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Risk = portfolio_risks,
    Return = portfolio_returns,
    Lambda = lambdas
  )
  
  # Calculate the return and risk of the specified portfolios
  portfolio_returns_list <- lapply(weights_list, function(weights) sum(mu * weights))
  portfolio_risks_list <- lapply(weights_list, function(weights) sqrt(t(weights) %*% Sigma %*% weights))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Risk, y = Return, color = Lambda)) +
    geom_line() +
    geom_point() +
    scale_color_gradient(low = "blue", high = "green") +
    labs(title = "Efficient Frontier",
         x = "Risk (Standard Deviation)",
         y = "Return",
         color = "Lambda") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # Add annotations for each portfolio
  portfolio_colors <- c("red", "dark blue", "dark green", "orange", "purple", "brown", "magenta", "black")
  for (i in 1:length(portfolios)) {
    p <- p + 
      annotate("point", x = portfolio_risks_list[[i]], y = portfolio_returns_list[[i]], 
               color = portfolio_colors[i], size = 3) +
      annotate("text", x = portfolio_risks_list[[i]], y = portfolio_returns_list[[i]],
               label = portfolios[i], vjust = -1.5, color = portfolio_colors[i])
  }
  
  # Print the plot
  print(p)
}

# Example usage
# Plot the entire dataset optimized portfolio weights
EWP <- EWP_portfolio_fun(my_dataset_list$`dataset 1`)
MVP <- MVP_portfolio_fun(my_dataset_list$`dataset 1`)
GMVP <- GMVP_portfolio_fun(my_dataset_list$`dataset 1`)
MSRP <- MSRP_portfolio_fun(my_dataset_list$`dataset 1`)
IVP <- IVP_portfolio_fun(my_dataset_list$`dataset 1`)
RPP <- RPP_portfolio_fun(my_dataset_list$`dataset 1`)
MDP <- MDP_portfolio_fun(my_dataset_list$`dataset 1`)
MDCP <- MDCP_portfolio_fun(my_dataset_list$`dataset 1`)
weights_list <- list(EWP,MVP,GMVP,MSRP,IVP,RPP,MDP,MDCP)
portfolios <- c("EWP","MVP","GMVP","MSRP","IVP","RPP","MDP","MDCP")
plot_efficient_frontier_many(my_dataset_list$`dataset 1`, weights_list, portfolios)

# Plot the backtest's rebalancing portfolio weights
EWP <- as.numeric(bt$EWP$`dataset 1`$w_optimized[3])
MVP <- as.numeric(bt$MVP$`dataset 1`$w_optimized[3])
GMVP <- as.numeric(bt$GMVP$`dataset 1`$w_optimized[3])
MSRP <- as.numeric(bt$MSRP$`dataset 1`$w_optimized[3])
IVP <- as.numeric(bt$IVP$`dataset 1`$w_optimized[3])
RPP <- as.numeric(bt$RPP$`dataset 1`$w_optimized[3])
MDP <- as.numeric(bt$MDP$`dataset 1`$w_optimized[3])
MDCP <- as.numeric(bt$MDCP$`dataset 1`$w_optimized[3])
weights_list <- list(EWP,MVP,GMVP,MSRP,IVP,RPP,MDP,MDCP)
portfolios <- c("EWP","MVP","GMVP","MSRP","IVP","RPP","MDP","MDCP")
plot_efficient_frontier_many(my_dataset_list$`dataset 1`, weights_list, portfolios)
