# Install and Load Required Packages
install.packages(c("tidyverse", "rstan", "gridExtra"))
library(tidyverse)
library(rstan)
library(gridExtra)

# Enable parallel processing for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load dataset
data1 <- read.csv("data_ecoli.csv")

# Scale observed infections to match N = 2200
case_S <- round(abs(jitter(data1$IS * (2200 / 187))))
case_R <- round(abs(jitter(data1$IR * (2200 / 187))))

# Total population
N <- 2200  

# Time points
n_days <- length(case_S)
t <- seq(1, n_days, by = 1)
t0 <- 0  

# Initial conditions
i_S0 <- 1
i_R0 <- 1
s0 <- N - i_S0 - i_R0
r0 <- 0
d1_0 <- 0
d2_0 <- 0
y0 <- c(S = s0, I_S = i_S0, I_R = i_R0, R = r0, D1 = d1_0, D2 = d2_0)

# Define the Stan Model
sir_model <- "
functions {
  real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real S = y[1];
    real I_S = y[2];
    real I_R = y[3];
    real R = y[4];
    real D1 = y[5];
    real D2 = y[6];
    real N = x_i[1];

    real lambda = theta[1];
    real beta = theta[2];
    real beta_c = theta[3];
    real delta_S = theta[4];
    real delta_R = theta[5];
    real mu = theta[6];
    real tau = theta[7];

    real dS_dt = lambda*N - beta * S * I_S / N - beta_c * S * I_R / N + mu * I_R + (mu + tau) * I_S;
    real dI_S_dt = beta * S * I_S / N - (mu + tau) * I_S - delta_S * I_S;
    real dI_R_dt = beta_c * S * I_R / N - mu * I_R - delta_R * I_R;
    real dR_dt = (mu + tau) * I_S + mu * I_R;
    real dD1_dt = delta_S * I_S;
    real dD2_dt = delta_R * I_R;

    return {dS_dt, dI_S_dt, dI_R_dt, dR_dt, dD1_dt, dD2_dt};
  }
}

data {
  int<lower=1> n_days;
  real y0[6];
  real t0;
  real ts[n_days];
  int N;
  int cases_S[n_days];
  int cases_R[n_days];
}

transformed data {
  real x_r[0];
  int x_i[1] = {N};
}

parameters {
  real<lower=0> lambda;
  real<lower=0> beta;
  real<lower=0> beta_c;
  real<lower=0> delta_S;
  real<lower=0> delta_R;
  real<lower=0> mu;
  real<lower=0> tau;
  real<lower=0> phi_inv;
}

transformed parameters {
  real y[n_days, 6];
  real phi = 1. / phi_inv;

  real theta[7];
  theta[1] = lambda;
  theta[2] = beta;
  theta[3] = beta_c;
  theta[4] = delta_S;
  theta[5] = delta_R;
  theta[6] = mu;
  theta[7] = tau;

  y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
}

model {
  // Priors
  lambda ~ exponential(1.544183e-08);
  beta ~ normal(0.06980519, 1);
  beta_c ~ normal(0.05769231, 1);
  delta_S ~ normal(0.04416404, 0.25);
  delta_R ~ normal(0.03597122, 0.25);
  mu ~ normal(0.03846154, 0.25);
  tau ~ normal(0.008075258, 0.25);
  phi_inv ~ exponential(5);

  // Likelihood
  cases_S ~ neg_binomial_2(y[, 2], phi);
  cases_R ~ neg_binomial_2(y[, 3], phi);
}

generated quantities {
  real recovery_time_S = 1 / (mu + tau);
  real recovery_time_R = 1 / mu;
  real pred_cases_S[n_days];
  real pred_cases_R[n_days];

  pred_cases_S = neg_binomial_2_rng(to_vector(y[, 2]) + rep_vector(1e-5, n_days), phi);
  pred_cases_R = neg_binomial_2_rng(to_vector(y[, 3]) + rep_vector(1e-5, n_days), phi);
}
"

# Prepare Data for Stan
data_sir <- list(
  n_days = n_days,
  y0 = y0,
  t0 = t0,
  ts = t,
  N = N,
  cases_S = case_S,
  cases_R = case_R
)

# Compile and Run the Stan Model
model <- stan_model(model_code = sir_model)

# Fit the model using MCMC sampling
niter <- 2000
fit_sir <- sampling(model, data = data_sir, iter = niter, chains = 3, seed = 0)

# Extract posterior parameters
pars <- c('lambda', 'beta', 'beta_c', 'delta_S', 'delta_R', "recovery_time_S", "recovery_time_R")
print(fit_sir, pars = pars)

# Summarize predictions
smr_pred <- cbind(
  as.data.frame(summary(fit_sir, pars = c("pred_cases_S", "pred_cases_R"), probs = c(0.05, 0.5, 0.95))$summary),
  t, case_S, case_R
)
colnames(smr_pred) <- make.names(colnames(smr_pred))  # Remove '%' from column names

smr_pred_S <- smr_pred[1:20,]
ggplot(smr_pred_S, aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "lightblue", alpha = 0.35) +  # Prediction uncertainty
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +  # Predicted I_S (median)
  geom_point(aes(y = case_S), color = "red", size = 1) +  # Observed I_S
  labs(x = "Day", y = "Number of Sensitive Infections (I_S)", 
       title = "Predicted vs Observed Sensitive Infections") +
  theme_minimal()




smr_pred_R <- smr_pred[21:40,]
ggplot(smr_pred_R, aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "lightblue", alpha = 0.35) +  # Prediction uncertainty
  geom_line(aes(y = X50.), color = "blue", linewidth = 1) +  # Predicted I_S (median)
  geom_point(aes(y = case_S), color = "red", size = 1) +  # Observed I_S
  labs(x = "Day", y = "Number of Sensitive Infections (I_S)", 
       title = "Predicted vs Observed Sensitive Infections") +
  theme_minimal()

