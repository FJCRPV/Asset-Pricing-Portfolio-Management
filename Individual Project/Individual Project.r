# -------------------------------------------------------------------------------------------
# Asset Pricing and Portfolio Management | INDIVIDUAL PROJECT | Francisco Perestrello, 20241560
# -------------------------------------------------------------------------------------------

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
p_load(dplyr)                 # for data manipulation
p_load(rcompanion)            # for statistical analysis
p_load(GLDEX)                 # for moments (mean, variance, skewness, kurtosis)
p_load(CVXR)                  # Convex Optimization in R
p_load(riskParityPortfolio)   # RPP
p_load(HierPortfolios)        # Hierarquical RPP
p_load(nloptr)                # Non-Convex Optimization in R
p_load(ggplot2)               # plot
p_load(portfolioBacktest)     # portfolio Backtesting

# -------------------------------------------------------------------------------------------
#                                       2. Individual Part
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# Stock price data parameter definitions
# -------------------------------------------------------------------------------------------

# Download data from Yahoo Finance
start_date <- "2010-01-01"
end_date   <- "2024-10-31"

Tickers <- c("SPY","GOOGL","NFLX","AAPL","AMZN","META","GLD","NVDA","IBM","TXN","ASML","DECK","V","MA","SNPS","JPM","AVGO","AMAT","TLT","NVO","LLY","AXON","XYL","GBTC","NTES","NEE","JNJ","PG","LMT","TSLA","LRCX")
# S&P 500, Bitcoin, Gold, 20+Yr US Bond and 27 Stocks ranging from Healthcare, to Financials, to Consumer electronic, Semiconductors, Software and Industrials

# -------------------------------------------------------------------------------------------
# Portfolio definitions
# -------------------------------------------------------------------------------------------

# Quintile portfolio
Quintile_portfolio_fun <- function(dataset, ...) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  N <- ncol(log_returns) 
  # design quintile portfolio
  ranking <- sort(colMeans(log_returns), decreasing = TRUE, index.return = TRUE)$ix
  w <- rep(0, N)
  w[ranking[1:round(N/5)]] <- 1/round(N/5)
  return(w)
}

# Hierarchical Risk Parity Portfolio (HRPP)                               
HRPP_portfolio_fun <- function(dataset, linkage = "single", graph = F, ...) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(log_returns)
  w <- HRP_Portfolio(Sigma, linkage = linkage, graph = graph) # The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
  return(as.vector(w$weights))
}

# Mean-Semivariance-Skewness Portfolio (MSVSP)
MSVSP_portfolio_fun <- function(dataset, gamma=1, lambda=0.5, ...) {
  # Step 1:  Calculate log returns, mean returns, and the variance-covariance matrix
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu    <- colMeans(log_returns)  # compute mean vector
  Sigma <- cov(log_returns)
  
  # Step 2: Calculate co-skewness terms
  t <- nrow(log_returns)
  n <- ncol(log_returns)
  coskewness <- array(0, dim = c(n, n, n))
  for(i in 1:n) {
    for(j in 1:n) {
      for(k in 1:n) {
          coskewness[i, j, k] <- sum(log_returns[,i] - mu[i] *
                                       log_returns[,j] - mu[j] *
                                       log_returns[,k] - mu[k]) / t
      }
    }
  }
  coskewness <- coskewness/sum(coskewness) # normalize
  coskewness_2d <- matrix(coskewness, nrow = n, ncol = n^2) # Convert to N * N^2 object
  
  # Step 3: Calculate co-semivariance terms
  cosemivariance <- matrix(0, nrow = n, ncol = n) 
  U <- which(apply(log_returns, 1, function(row) any(row < 0))) # Indices where returns were negative
  for (i in 1:n) {
    for (j in 1:n) {
      cosemivariance[i, j] <- sum((log_returns[U, i]) * (log_returns[U, j])) / t
    }
  }
  cosemivariance <- cosemivariance/sum(cosemivariance) # normalize
  
  # Step 3: Set up the objective function
  objective_function <- function(w) {
    # Calculate portfolio skewness
    portfolio_skewness <- t(w) %*% coskewness_2d %*% kronecker(w, w)
    
    # Calculate portfolio semivariance
    portfolio_semivariance <- t(w) %*% cosemivariance %*% w
    
    # Objective function
    objective_value <- gamma * portfolio_skewness - lambda * portfolio_semivariance
    
    # Maximize by returning the negative for minimization solver
    return(-objective_value)
  }
  
  # Constraints
  constraints_function <- function(w) {
    return(sum(w) - 1)
  }
  
  # Initial values
  w0 <- rep(1 / n, n)  # Equal weights as starting point
  
  # Optimization with nloptr
  result <- nloptr(
    x0 = w0,
    eval_f = objective_function,
    lb = rep(0, n),  # Long-only constraint
    ub = rep(1, n),  # Upper bound for weights
    eval_g_eq = constraints_function,
    opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 10000)         # Note: this algorithm uses Linear Approximations to resolve a non-convex problem
  )                                                                     # https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  
  # Get optimized weights
  return(result$solution)
}
  
# -------------------------------------------------------------------------------------------
# Backtesting
# -------------------------------------------------------------------------------------------

prices <- stockDataDownload(stock_symbols = Tickers, 
                            index_symbol = "^GSPC",
                            from = start_date, to = end_date)

save(prices, file = "stockdata_from_2010-01-01_to_today_individual_part.Rdata")
#load('stockdata_from_2010-01-01_to_today_individual_part.Rdata')

# Resample 100 times, each with 10 stocks and 2-year consecutive data
set.seed(123)

my_dataset_list <- financialDataResample(prices, 
                                         N_sample = 10,     # Desired number of financial instruments in each resample
                                         T_sample = 252*2,  # Desired length of each resample 
                                         num_datasets = 100) # Number of resampled datasets

portfolios <- list("Quintile" = Quintile_portfolio_fun,
                   "HRPP" = HRPP_portfolio_fun,
                   "MSVSP" = MSVSP_portfolio_fun)

bt <- portfolioBacktest(portfolios, my_dataset_list,
                        benchmark = ("index"),
                        rebalance_every = 252*0.5,           # arbitrary
                        optimize_every = 252*0.5,            # arbitrary
                        lookback = 252*0.5,                  # arbitrary
                        shortselling = F,                 
                        cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4))  # arbitrary

# Summary of portfolios performance measures
res_sum <- backtestSummary(bt, summary_fun = median, show_benchmark = TRUE) # Show median values on the table

summaryTable(res_sum, type = "DT", 
             order_col = "Sharpe ratio",   # performance measure to be used to sort the rows
             order_dir = "desc")

# BoxPlot of portfolio simulations' Sharpe Ratios:
backtestBoxPlot(bt, measure = "Sharpe ratio")

# Cumulative Return or wealth plot of a single backtest
backtestChartCumReturn(bt, dataset_num = 1)

# Drawdown plot of a single backtest
backtestChartDrawdown(bt, dataset_num = 1)

# Portfolio allocation evolution of a particular portfolio over a single backtest
backtestChartStackedBar(bt, "MSVSP", dataset_num = 1, legend = TRUE)

# -------------------------------------------------------------------------------------------
# Strategy Analysis
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# Hierarchical Risk Parity Portfolio (HRPP)
# -------------------------------------------------------------------------------------------

# Plot the dendrogram
weights <- HRPP_portfolio_fun(my_dataset_list$`dataset 1`, graph = T)

# Function to plot correlation matrix heatmap
plot_correlation_heatmap <- function(dataset) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cor(log_returns)  # compute correlation matrix  
  
  # Convert correlation matrix to data frame for ggplot
  cov_df <- as.data.frame(as.table(Sigma))
  names(cov_df) <- c("Asset1", "Asset2", "Correlation")
  
  # Create heatmap using ggplot2
  p <- ggplot(cov_df, aes(x = Asset1, y = Asset2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         midpoint = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle("Correlation Matrix Heatmap")
  
  return(p)
}
# Example usage
plot_correlation_heatmap(my_dataset_list$`dataset 1`)


# Function to plot clustered correlation matrix
plot_clustered_correlation_heatmap <- function(dataset) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cor(log_returns)  # compute correlation matrix
  clust <- hclust(dist(Sigma), method = 'single')
  clustOrder <- clust$order
  
  # Re-order the correlation matrix
  clustured_Sigma <- Sigma[clustOrder, clustOrder]
  
  # Convert correlation matrix to data frame for ggplot
  cluster_cor_df <- as.data.frame(as.table(clustured_Sigma))
  names(cluster_cor_df) <- c("Asset1", "Asset2", "Correlation")
  
  # Create heatmap using ggplot2
  p <- ggplot(cluster_cor_df, aes(x = Asset1, y = Asset2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         midpoint = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle("Clustered Correlation Matrix Heatmap")
  
  return(p)
}

# Example usage
plot_clustered_correlation_heatmap(my_dataset_list$`dataset 1`)


# Function to plot covariance matrix heatmap
plot_covariance_heatmap <- function(dataset) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(log_returns)  # compute covariance matrix  
  
  # Convert covariance matrix to data frame for ggplot
  cov_df <- as.data.frame(as.table(Sigma))
  names(cov_df) <- c("Asset1", "Asset2", "Covariance")
  
  # Create heatmap using ggplot2
  p <- ggplot(cov_df, aes(x = Asset1, y = Asset2, fill = Covariance)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         midpoint = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle("Covariance Matrix Heatmap")
  
  return(p)
}

# Example usage
plot_covariance_heatmap(my_dataset_list$`dataset 1`)


# Function to plot clustered covariance matrix
plot_clustered_covariance_heatmap <- function(dataset) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  corr <- cor(log_returns)  # compute correlation matrix
  clust <- hclust(dist(corr), method = 'single')
  clustOrder <- clust$order
  
  # Re-order the correlation matrix
  Sigma <- cov(log_returns)
  clustured_Sigma <- Sigma[clustOrder, clustOrder]
  
  # Convert correlation matrix to data frame for ggplot
  cluster_cov_df <- as.data.frame(as.table(clustured_Sigma))
  names(cluster_cov_df) <- c("Asset1", "Asset2", "Covariance")
  
  # Create heatmap using ggplot2
  p <- ggplot(cluster_cov_df, aes(x = Asset1, y = Asset2, fill = Covariance)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         midpoint = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle("Clustered Covariance Matrix Heatmap")
  
  return(p)
}

# Example usage
plot_clustered_covariance_heatmap(my_dataset_list$`dataset 1`)


# Function to plot weights and relative risk contributions
plot_w_rrc <- function(dataset, weights) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(log_returns)  # compute covariance matrix 
  
  # Plot
  barplotPortfolioRisk(weights, Sigma) +
    ggtitle('Hierarchical Risk Parity Portfolio (HRPP)') +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = colnames(Sigma))
}

# Example usage
weights <- HRPP_portfolio_fun(my_dataset_list$`dataset 1`)
plot_w_rrc(my_dataset_list$`dataset 1`, weights)

# -------------------------------------------------------------------------------------------
# Mean-Semivariance Skewness Portfolio (MSVSP)
# -------------------------------------------------------------------------------------------

# Skewness Semivariance Efficient Frontier
plot_skew_semivar_efficient_frontier <- function(dataset, weights, portfolio) {
  log_returns <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu <- colMeans(log_returns)  # compute mean vector
  Sigma <- cov(log_returns)  # compute covariance matrix
  
  # Calculate co-skewness and co-semivariance terms
  t <- nrow(log_returns)
  n <- ncol(log_returns)
  coskewness <- array(0, dim = c(n, n, n))
  for(i in 1:n) {
    for(j in 1:n) {
      for(k in 1:n) {
        coskewness[i, j, k] <- sum(log_returns[,i] - mu[i] *
                                     log_returns[,j] - mu[j] *
                                     log_returns[,k] - mu[k]) / t
      }
    }
  }
  coskewness <- coskewness/sum(coskewness) # normalize
  coskewness_2d <- matrix(coskewness, nrow = n, ncol = n^2)
  
  cosemivariance <- matrix(0, nrow = n, ncol = n) 
  U <- which(apply(log_returns, 1, function(row) any(row < 0))) # Indices where returns were negative
  for (i in 1:n) {
    for (j in 1:n) {
      cosemivariance[i, j] <- sum((log_returns[U, i]) * (log_returns[U, j])) / t
    }
  }
  cosemivariance <- cosemivariance/sum(cosemivariance) # normalize
  
  # Create a range of risk aversion parameters
  lambdas <- seq(0, 1, length.out = 100)
  
  # Initialize vectors to store portfolio metrics
  portfolio_skewnesses <- numeric(length(lambdas))
  portfolio_semivariances <- numeric(length(lambdas))
  
  # Calculate portfolios for different risk aversion parameters
  for(i in seq_along(lambdas)) {
    # Get the optimized weights
    w_opt <- MSVSP_portfolio_fun(dataset, lambda = lambdas[i]) # Use the skewness-semivariance portfolio to draw the frontier
    
    # Compute portfolio metrics
    portfolio_skewnesses[i] <- t(w_opt) %*% coskewness_2d %*% kronecker(w_opt, w_opt)
    portfolio_semivariances[i] <- t(w_opt) %*% cosemivariance %*% w_opt
  }
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Semivariance = portfolio_semivariances,
    Skewness = portfolio_skewnesses,
    Lambda = lambdas
  )

  # Calculate the skewness and semivariance of the specified portfolio
  portfolio_skewness <- t(weights) %*% coskewness_2d %*% kronecker(weights, weights)
  portfolio_semivariance <- t(weights) %*% cosemivariance %*% weights
  
  # Create the plot
  ggplot(plot_data, aes(x = Semivariance, y = Skewness, color = Lambda)) +
    geom_line() +
    geom_point() +
    scale_color_gradient(low = "blue", high = "green") +
    labs(title = "Skewness-Semivariance Efficient Frontier",
         x = "Semivariance",
         y = "Skewness",
         color = "Lambda") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    ) +
    annotate("point", x = portfolio_semivariance, y = portfolio_skewness, 
             color = "black", size = 3) +
    annotate("text", x = portfolio_semivariance, y = portfolio_skewness,
             label = portfolio, vjust = -1.5, color = "black")
}

# Example usage:
weights <- MSVSP_portfolio_fun(my_dataset_list$`dataset 1`)
plot_skew_semivar_efficient_frontier(my_dataset_list$`dataset 1`, weights, "MSVSP")
# To plot a zoomed-out version, one should change the values in the lambdas sequence inside the function

# Plot the backtest's optimized portfolio weights on the frontier
opt_weights <- as.numeric(bt$MSVSP$`dataset 1`$w_optimized[3])
plot_skew_semivar_efficient_frontier(my_dataset_list$`dataset 1`, opt_weights, "MSVSP")


# Function to plot asset skewness vs semivariance, with point width dependent on the weights
plot_skewness_vs_semivariance <- function(dataset, weights, title = "Skewness vs. Semivariance") {
  log_returns <- diff(log(dataset$adjusted))[-1]
  mu <- colMeans(log_returns)
  
  # Calculate skewnesses and variances
  semi_variances <- colMeans(pmin(log_returns - matrix(rep(mu, nrow(log_returns)), nrow=nrow(log_returns), byrow=TRUE), 0))
  skewness_values <- apply(log_returns, 2, skewness)
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Asset = colnames(log_returns),
    Semivariance = semi_variances,
    Skewness = skewness_values,
    Weight = weights
  )
  
  # Create the bubble chart
  ggplot(plot_data, aes(x = Semivariance, y = Skewness, size = Weight, label = Asset)) +
    geom_point(alpha = 0.5) +
    geom_text(vjust = 1.5, hjust = 0.5, check_overlap = TRUE) +  # Adjusted position to ensure all labels are visible
    labs(title = title, x = "Semivariance", y = "Skewness") +
    theme_minimal() +
    scale_size_continuous(range = c(1, 10)) +  # Adjust bubble size range
    xlim(min(semi_variances) - 0.001, max(semi_variances) + 0.001) +  # Ensure all stocks are present in the plot
    ylim(min(skewness_values) - 0.01, max(skewness_values) + 0.01)  # Ensure all stocks are present in the plot
}

# Example usage
weights <- MSVSP_portfolio_fun(my_dataset_list$`dataset 1`) # Weights com a função
plot_skewness_vs_semivariance(my_dataset_list$`dataset 1`, weights)

# -------------------------------------------------------------------------------------------
# Mean-Variance Efficient Frontier Analysis
# -------------------------------------------------------------------------------------------

# Mean-Variance Portfolio
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

# Function to plot strategies on the Markowitz Mean-Variance Efficient Frontier
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
    w_opt <- MVP_portfolio_fun(dataset, lambda = lambdas[i])
    
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
  portfolio_colors <- c("red", "orange", "purple")
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
Quintile <- Quintile_portfolio_fun(my_dataset_list$`dataset 1`)
HRPP <- HRPP_portfolio_fun(my_dataset_list$`dataset 1`)
MSVSP <- MSVSP_portfolio_fun(my_dataset_list$`dataset 1`)
weights_list <- list(Quintile,HRPP,MSVSP)
portfolios <- c("Quintile","HRPP","MSVSP")
plot_efficient_frontier_many(my_dataset_list$`dataset 1`, weights_list, portfolios)

# Plot the backtest optimized portfolio weights
Quintile <- as.numeric(bt$Quintile$`dataset 1`$w_optimized[3])
HRPP <- as.numeric(bt$HRPP$`dataset 1`$w_optimized[3])
MSVSP <- as.numeric(bt$MSVSP$`dataset 1`$w_optimized[3])
weights_list <- list(Quintile,HRPP,MSVSP)
portfolios <- c("Quintile","HRPP","MSVSP")
plot_efficient_frontier_many(my_dataset_list$`dataset 1`, weights_list, portfolios)

