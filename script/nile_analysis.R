# Nile Time Series Analysis - Local Level Model
# Authors: Andrea Bianco, Egle Caudullo, Giovanni Contento
# Date: May 2025

# 1. Libraries
library(TSA)
library(dlm)
library(statespacer)
library(graphics)
library(stats)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(float)
library(tidyverse)

# 2. Data Loading
data(Nile)
y <- as.numeric(Nile)
n <- length(y)
Nile_mat <- matrix(Nile, ncol = 1)
Nile_df <- data.frame(year = time(Nile), flow = as.numeric(Nile))

# Descriptive statistics
cat("Descriptive Statistics:\n")
summary(Nile)
cat("\nStandard Deviation:", sd(Nile), "\n")

# 3. Data Visualization
par(mfrow = c(1, 1))
plot(Nile, 
     main = "Nile time series", 
     xlab = "Time", 
     ylab = "Annual flow volume",
     col = "darkorange",
     lwd = 2)

# 4. Kalman Filter - Model Fitting
fit <- statespacer(Nile_mat,
                   local_level_ind = TRUE,
                   diagnostics = TRUE,
                   initial = 0.5 * log(var(Nile_mat)))

# Extract variance parameters
sigma2_eps <- as.vector(fit$system_matrices$H$H)
sigma2_eta <- as.vector(fit$system_matrices$Q$full)

cat("\n=== Estimated Parameters ===\n")
cat("Observation variance (σ²_ε):", sigma2_eps, "\n")
cat("State variance (σ²_η):", sigma2_eta, "\n")
cat("Signal-to-noise ratio (q):", sigma2_eta / sigma2_eps, "\n")

# 5. Kalman Filter Output Visualization
filtered_state <- fit$filtered$level
filtered_variance <- fit$filtered$P
prediction_errors <- fit$diagnostics$v
prediction_variance <- fit$diagnostics$F

par(mfrow = c(2, 2))

# Panel 1: Filtered level with confidence intervals
plot(time(Nile), filtered_state, type = "l", lwd = 2,
     main = "Filtered level with 90% CI, and observed data points",
     xlab = "", ylab = "")
points(time(Nile), Nile, pch = 19, cex = 0.5)
ci_width <- 1.645 * sqrt(filtered_variance)
lines(time(Nile), filtered_state + ci_width, lty = 2, col = "red")
lines(time(Nile), filtered_state - ci_width, lty = 2, col = "red")

# Panel 2: Filtered state variance
plot(time(Nile), filtered_variance, type = "l", lwd = 2,
     main = "Filtered state variance P_t",
     xlab = "", ylab = "")

# Panel 3: Prediction errors
plot(time(Nile), prediction_errors, type = "l", lwd = 2,
     main = "Predictions errors v_t",
     xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "gray")

# Panel 4: Prediction variance
plot(time(Nile), prediction_variance, type = "l", lwd = 2,
     main = "Predictions variance F_t",
     xlab = "", ylab = "")

par(mfrow = c(1, 1))

# 6. Kalman Smoothing
smoothed_state <- fit$smoothed$level
smoothed_variance <- fit$smoothed$V

# Compute smoothed disturbances with standard deviations
D_t <- as.vector(fit$diagnostics$D)
N_t <- as.vector(fit$diagnostics$N)

# Observation disturbance
smoothed_obs_dist <- fit$smoothed$epsilon
var_obs_dist <- sigma2_eps - (sigma2_eps^2) * D_t
std_obs_dist <- sqrt(var_obs_dist)

# State disturbance
smoothed_state_dist <- fit$smoothed$eta
var_state_dist <- sigma2_eta - (sigma2_eta^2) * N_t
std_state_dist <- sqrt(var_state_dist)

# 7. Smoothing Output Visualization
par(mfrow = c(2, 2))

# Smoothed observation disturbances
plot(time(Nile), smoothed_obs_dist, type = "l", lwd = 2,
     main = "Smoothed Observation Disturbances",
     xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "gray")

# Observation error standard deviation
plot(time(Nile), std_obs_dist, type = "l", lwd = 2,
     main = "Observation Error Standard Deviation",
     xlab = "", ylab = "")

# Smoothed state disturbances
plot(time(Nile), smoothed_state_dist, type = "l", lwd = 2,
     main = "Smoothed State Disturbances",
     xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "gray")

# State error standard deviation
plot(time(Nile), std_state_dist, type = "l", lwd = 2,
     main = "State Error Standard Deviation",
     xlab = "", ylab = "")

par(mfrow = c(1, 1))

# 8. Comparison: Filtered vs Smoothed
plot(time(Nile), Nile, type = "l", lwd = 1.5, col = "black",
     main = "", xlab = "Time", ylab = "",
     ylim = range(c(Nile, filtered_state, smoothed_state)))
lines(time(Nile), filtered_state, col = "red", lwd = 2)
lines(time(Nile), smoothed_state, col = "blue", lwd = 2)
legend("topright", 
       legend = c("Observed", "Filtered", "Smoothed"),
       col = c("black", "red", "blue"), 
       lwd = c(1.5, 2, 2))

# 9. Simulation
a1 <- Nile[1]
set.seed(57)

# Generate random disturbances
eta_t_sim <- rnorm(n, 0, sqrt(sigma2_eta))
epsilon_t_sim <- rnorm(n, 0, sqrt(sigma2_eps))

# Simulate observations
y_sim <- a1 + cumsum(eta_t_sim) + epsilon_t_sim
alpha_t_sim <- y_sim - epsilon_t_sim

# Fit model to simulated data
y_mat <- matrix(y_sim)
fit_sim <- statespacer(y_mat,
                       local_level_ind = TRUE,
                       initial = 0.5 * log(var(y_mat)))

# Conditional simulation
epsilon_tilde <- epsilon_t_sim - fit_sim$smoothed$epsilon + fit$smoothed$epsilon
alpha_tilde <- Nile_df$flow - epsilon_tilde
alpha_tilde_1 <- as.vector(alpha_tilde[2:length(alpha_tilde)])
alpha_tilde_1[100] <- 774.9514
eta_tilde <- alpha_tilde_1 - alpha_tilde

# 10. Simulation Visualization
par(mfrow = c(2, 2))

# Unconditional sample
plot(time(Nile), alpha_t_sim, pch = 1, cex = 0.8,
     main = "sample state (dots)\nsmoothed state (solid line)",
     xlab = "", ylab = "")
lines(time(Nile), smoothed_state, lwd = 2)

# Conditional sample state
plot(time(Nile), alpha_tilde, pch = 19, cex = 0.6,
     main = "conditional sample state(dots)\nsmoothed state (solid line)",
     xlab = "", ylab = "")
lines(time(Nile), smoothed_state, lwd = 2)

# Conditional sample epsilon
plot(time(Nile), epsilon_tilde, pch = 19, cex = 0.6,
     main = "conditional sample epsilon (dots)\nsmoothed observation error (solid line)",
     xlab = "", ylab = "")
lines(time(Nile), smoothed_obs_dist, lwd = 2)

# Conditional sample eta
plot(time(Nile), eta_tilde, pch = 19, cex = 0.6,
     main = "conditional sample eta(dots)\nsmoothed state error (solid line)",
     xlab = "", ylab = "")
lines(time(Nile), smoothed_state_dist, lwd = 2)

par(mfrow = c(1, 1))

# 11. Parameter Estimation - Concentrated Diffuse Log-Likelihood
CDLL <- function(y, psi) {
  a2 <- y[1]
  P2 <- 1 + exp(psi)
  n <- length(y)
  
  # Initialize vectors
  a_pred <- numeric(n)
  P <- numeric(n)
  v <- numeric(n)
  K <- numeric(n)
  F <- numeric(n)
  dllk <- numeric(n)
  
  a_pred[2] <- a2
  P[2] <- P2
  
  # Kalman filter recursion
  for (t in 2:(n - 1)) {
    v[t] <- y[t] - a_pred[t]
    F[t] <- P[t] + 1
    K[t] <- P[t] / F[t]
    P[t + 1] <- P[t] * (1 - K[t]) + exp(psi)
    a_pred[t + 1] <- a_pred[t] + K[t] * v[t]
    dllk[t] <- -0.5 * log(F[t])
  }
  
  # Final observation
  F[n] <- P[n] + 1
  K[n] <- P[n] / F[n]
  dllk[n] <- -0.5 * log(F[n])
  v[n] <- y[n] - a_pred[n]
  
  # Concentrated estimate
  sigma_e_hat <- sum(v[2:n]^2 / F[2:n]) / (n - 1)
  llk <- -0.5 * (n - 1) * log(sigma_e_hat) + sum(dllk[2:n])
  
  return(list(a_pred = a_pred, 
              llk = llk, 
              v = v, 
              F = F, 
              sigma_e_hat = sigma_e_hat))
}

# Log-likelihood function with recording
loglist <- list()

loglikelihood <- function(par, y, record = TRUE) {
  psi <- par[1]
  
  if (exp(psi) <= 0) return(1e10)
  
  dll_out <- CDLL(y, psi)
  llk <- dll_out$llk
  
  if (record) {
    h <- 1e-6
    l_plus <- loglikelihood(psi + h, y, record = FALSE)
    l_minus <- loglikelihood(psi, y, record = FALSE)
    score <- (l_plus - l_minus) / h
    
    loglist[[length(loglist) + 1]] <<- list(
      Iteration = length(loglist),
      psi = psi,
      logLik = llk,
      score = score
    )
  }
  
  return(-llk)
}

# Optimizer function
optimizer <- function(y, par, method = "BFGS") {
  loglist <<- list()
  optim(par = par, 
        fn = loglikelihood, 
        y = y, 
        method = method,
        control = list(reltol = 1e-12, maxit = 50))
}

# 12. Run Optimization
q_0 <- 1
psi_0 <- log(q_0)

# Run optimization
res <- optimizer(y, psi_0, method = "BFGS")

# Extract results
results_df <- data.frame(
  Iteration = sapply(loglist, function(x) x$Iteration),
  q = sapply(loglist, function(x) exp(round(x$psi, 4))),
  psi = sapply(loglist, function(x) round(x$psi, 2)),
  score = sapply(loglist, function(x) round(x$score, 5)),
  logLikelihood = sapply(loglist, function(x) round(x$logLik, 2))
)

cat("\n=== Optimization Results ===\n")
print(results_df)

cat("\n=== Final Estimates ===\n")
cat("Optimal q:", exp(res$par), "\n")
cat("Optimal ψ:", res$par, "\n")
cat("Log-Likelihood:", -res$value, "\n")
cat("Convergence:", res$convergence, "(0 = success)\n")
