livingexp.dat<-read.csv("/Users/tomasmock/LSE/Bayesian Data Analysis/Reading week project/LES_RW.csv",header = TRUE, stringsAsFactors = TRUE)

library(gridExtra)
library(ggplot2)
library(tidyverse)
library(arm)
library(nimble)
library(rstan)
library(rstanarm)
library(bayesplot)
library(loo) 
library(dplyr)

# Data preprocessing #

livingexp.dat <- dplyr::select(livingexp.dat, c("expenditure","income","Gorx","A121r"))%>% mutate(A121r=factor(A121r))%>%
  mutate(income=ifelse(income==0,0.5,income))%>%
  mutate(expenditure=log(expenditure)) %>%
  mutate(income=log(income))

no_low_livingexp.dat <- livingexp.dat %>% filter(income>2.5)

# centre the age and income variables
livingexp.dat.tr <- no_low_livingexp.dat %>%
  mutate(expenditure = expenditure - mean(expenditure)) %>%
  mutate(income = income - mean(income))


# Model 1 (Random slope + Intercept with STAN) #

# Read and prepare your data
livingexp.dat.tr <- livingexp.dat.tr %>%
  mutate(A121r = relevel(A121r, "1"))


# Create model matrix for fixed effects
X <- model.matrix(expenditure ~ income + A121r, data = livingexp.dat.tr)

print(colnames(X))


# Ensure that Gorx is a factor and convert to numeric
livingexp.dat.tr$Gorx <- as.numeric(as.factor(livingexp.dat.tr$Gorx))

# Define the number of predictors and hierarchical levels
D <- ncol(X)
nGorx <- length(unique(livingexp.dat.tr$Gorx))

# Prepare the data list for STAN
hierarchical_stan_data <- list(
  expenditure = livingexp.dat.tr$expenditure,
  X = X,
  N = nrow(livingexp.dat.tr),
  D = D,
  Gorx = livingexp.dat.tr$Gorx,
  nGorx = nGorx
)

# Define initial values for two chains
hierarchical_stan_inits <- list(
  list(beta = rep(0, D), nu_intercept = rep(0, nGorx), nu_slope = rep(0, nGorx), sigma_res = 1, sigma_nu_intercept = 1, sigma_nu_slope = 1),
  list(beta = rep(1, D), nu_intercept = rep(1, nGorx), nu_slope = rep(1, nGorx), sigma_res = 10, sigma_nu_intercept = 10, sigma_nu_slope = 10)
)

# Write the STAN model to a file
writeLines(con = "stan_models/hierarchical_living.stan", text = "
data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> nGorx;
  matrix[N, D] X;
  int<lower=1> Gorx[N];
  vector[N] expenditure;
}

parameters {
  vector[D] beta;
  real<lower=0> sigma_res;
  vector[nGorx] nu_intercept;
  vector[nGorx] nu_slope;
  real<lower=0> sigma_nu_intercept;
  real<lower=0> sigma_nu_slope;
}

model {
  beta ~ normal(0, 2.5);
  nu_intercept ~ normal(0, sigma_nu_intercept);
  nu_slope ~ normal(0, sigma_nu_slope);
  sigma_res ~ normal(0, 5);
  sigma_nu_intercept ~ normal(0, 5);
  sigma_nu_slope ~ normal(0, 5);
  for (i in 1:N) {
    expenditure[i] ~ normal(dot_product(beta, X[i]) + nu_intercept[Gorx[i]] + nu_slope[Gorx[i]] * X[i, 2], sigma_res);
  }
}

generated quantities {
  vector[N] expenditure_pred;
  vector[N] log_lik;
  for (i in 1:N) {
    expenditure_pred[i] = normal_rng(dot_product(beta, X[i]) + nu_intercept[Gorx[i]] + nu_slope[Gorx[i]] * X[i, 2], sigma_res);
    log_lik[i] = normal_lpdf(expenditure[i] | dot_product(beta, X[i]) + nu_intercept[Gorx[i]] + nu_slope[Gorx[i]] * X[i, 2], sigma_res);
  }
}
")

# Compile and sample from the model
hierarchical_stan_model <- stan_model("stan_models/hierarchical_living.stan", model_name = "hierarchical_living")
hierarchical_living_stan_fit <- sampling(hierarchical_stan_model, data = hierarchical_stan_data, warmup = 4000, iter = 10000, chains = 2, thin = 2, init = hierarchical_stan_inits)



### summary as an object ###
fit_summary <- summary(hierarchical_living_stan_fit)

# view the summary
print(fit_summary)




# Extract predictions
expenditure_pred_extracted <- extract(hierarchical_living_stan_fit)$expenditure_pred

# visualize the entire distribution of predictions
ppc_dens_overlay(y = livingexp.dat.tr$expenditure, yrep = expenditure_pred_extracted)





# Calculate the summary statistics for observed and predicted expenditures
# Here we calculate both mean and standard deviation
obs_mean <- mean(livingexp.dat.tr$expenditure)
obs_sd <- sd(livingexp.dat.tr$expenditure)

pred_mean <- apply(expenditure_pred_extracted, 2, mean)
pred_sd <- apply(expenditure_pred_extracted, 2, sd)

# Using ppc_stat_2d to compare the joint distributions of mean and sd of observed and predicted expenditures
ppc_stat_2d(
  y = livingexp.dat.tr$expenditure,
  yrep = expenditure_pred_extracted,
  stat = c("mean", "sd")
)






# Extract the log likelihood matrix from the Stan fit
log_lik <- extract_log_lik(hierarchical_living_stan_fit, parameter_name = "log_lik", merge_chains = FALSE)

# Compute the relative efficiency of the importance sampling estimator
r_eff <- relative_eff(exp(log_lik), cores = 2)

# Compute the LOO using the loo package
livingexp_loo <- loo(log_lik, r_eff = r_eff, cores = 2, save_psis = TRUE)

# Print the LOO output
print(livingexp_loo)




# WAIC #
# Calculate WAIC 
livingexp_waic <- waic(log_lik, r_eff = r_eff, cores = 2)

# Print the WAIC output
print(livingexp_waic)




# LOO-PIT #
psis <- livingexp_loo$psis_object
lw <- weights(psis)


#  'expenditure_pred' extracted from hierarchical Stan model
expenditure_pred <- extract(hierarchical_living_stan_fit, pars = "expenditure_pred")$expenditure_pred

# LOO-PIT Overlay plot
ppc_loo_pit_overlay(y = livingexp.dat.tr$expenditure, yrep = expenditure_pred, lw = lw)






# Extract summary statistics for random intercepts
summary_extract <- summary(hierarchical_living_stan_fit, pars = paste("nu_intercept[", 1:12, "]", sep = ""))
ri_cdcnt <- as.data.frame(summary_extract$summary)[, c(1, 4, 8)]
colnames(ri_cdcnt) <- c("mean", "lowci", "upci")

# Prepare the data frame for plotting
forplot_cdcnts <- data.frame(
  region = paste("region", 1:12, sep = "."),
  mean = ri_cdcnt$mean,
  lowci = ri_cdcnt$lowci,
  upci = ri_cdcnt$upci
)

# Arrange the data by mean for better visualization
forplot_cdcnts <- arrange(forplot_cdcnts, mean)

# Add a rank column to order the patients in the plot
forplot_cdcnts <- mutate(forplot_cdcnts, pat_rank = seq(1, nrow(forplot_cdcnts)))

# Create the forest plot using ggplot2
fp <- ggplot(data = forplot_cdcnts, aes(x = region, y = mean, ymin = lowci, ymax = upci)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  labs(x = "Region", y = "Intercept Estimate") +
  coord_flip()  # Flips the axes for easier reading

# Display the plot
print(fp)




# Extract summary statistics for random slope
summary_extract <- summary(hierarchical_living_stan_fit, pars = paste("nu_slope[", 1:12, "]", sep = ""))
ri_cdcnt <- as.data.frame(summary_extract$summary)[, c(1, 4, 8)]
colnames(ri_cdcnt) <- c("mean", "lowci", "upci")

# Prepare the data frame for plotting
forplot_cdcnts <- data.frame(
  region = paste("region", 1:12, sep = "."),
  mean = ri_cdcnt$mean,
  lowci = ri_cdcnt$lowci,
  upci = ri_cdcnt$upci
)

# Arrange the data by mean for better visualization
forplot_cdcnts <- arrange(forplot_cdcnts, mean)

# Add a rank column to order the patients in the plot
forplot_cdcnts <- mutate(forplot_cdcnts, pat_rank = seq(1, nrow(forplot_cdcnts)))

# Create the forest plot using ggplot2
fp <- ggplot(data = forplot_cdcnts, aes(x = region, y = mean, ymin = lowci, ymax = upci)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  labs(x = "Region", y = "Slope Estimate") +
  coord_flip()  # Flips the axes for easier reading

# Display the plot
print(fp)



# More plots #

bayesplot::mcmc_pairs(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]'))


bayesplot::mcmc_acf(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]'))


bayesplot::mcmc_trace(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]'))







# MODEL 1 Hierarchical modelling with NIMBLE (Random Intercept + Slope)

# Gorx as a factor
livingexp.dat.tr$Gorx <- as.factor(livingexp.dat.tr$Gorx)

summary(livingexp.dat.tr)


# Model setup
nGroups <- length(unique(livingexp.dat.tr$Gorx)) # Number of groups
N <- nrow(livingexp.dat.tr)                      # Number of observations
J <- length(unique(livingexp.dat.tr$A121r))      # Number of levels in A121r

# Create a NIMBLE model
nim_int_slp <- nimbleCode({
  beta0 ~ dnorm(0, 0.25) 
  betaIncome ~ dnorm(0, 0.25)
  for (j in 1:J) {
    betaA121r[j] ~ dnorm(0, 0.25)
  }
  for (g in 1:nGroups) {
    interceptGorx[g] ~ dnorm(0, tau_gorx)
    slopeGorx[g] ~ dnorm(0, tau_gorx)
  }
  tau_gorx ~ dgamma(2.5, 0.5)
  sigma_gorx <- sqrt(1 / tau_gorx)
  tau_res ~ dgamma(2.5, 0.5)
  sigma_res <- sqrt(1 / tau_res)
  for (i in 1:N) {
    mu[i] <- beta0 + betaIncome * income[i] +
      inprod(betaA121r[1:J], A121r[i,1:J]) +
      interceptGorx[GorxID[i]] +
      slopeGorx[GorxID[i]] * income[i]
    expenditure[i] ~ dnorm(mu[i], tau_res)
  }
})

# Prepare data and initial values for NIMBLE
A121r_matrix <- model.matrix(~ A121r - 1, data = livingexp.dat.tr)
Gorx_indices <- as.numeric(livingexp.dat.tr$Gorx)

nimble_data <- list(expenditure = livingexp.dat.tr$expenditure,
                    income = livingexp.dat.tr$income,
                    A121r = A121r_matrix,
                    GorxID = Gorx_indices)

nimble_constants <- list(N = N, J = J, nGroups = nGroups)
nimble_inits <- list(beta0 = 0, betaIncome = 0, betaA121r = rep(0, J),
                     interceptGorx = rep(0, nGroups), slopeGorx = rep(0, nGroups),
                     tau_gorx = 1, tau_res = 1)

# Build and compile the model, enabling WAIC
nim_int_slp_model <- nimbleModel(code = nim_int_slp, 
                                 constants = nimble_constants, 
                                 data = nimble_data, 
                                 inits = nimble_inits)
compiled_nim_int_slp_model <- compileNimble(nim_int_slp_model)
nim_int_slp_mcmc <- buildMCMC(nim_int_slp_model, enableWAIC = TRUE)
compiled_nim_int_slp_mcmc <- compileNimble(nim_int_slp_mcmc, project = nim_int_slp_model)

# Run the model
m <- 10000
bi <- 4000
nimble_inits <- list(
  list(beta0 = 0, betaIncome = 0.1, betaA121r = rep(0.1, J), interceptGorx = rep(0.1, nGroups),
       slopeGorx = rep(0.1, nGroups), tau_gorx = 1, tau_res = 1),
  list(beta0 = 0, betaIncome = -0.1, betaA121r = rep(-0.1, J), interceptGorx = rep(-0.1, nGroups),
       slopeGorx = rep(-0.1, nGroups), tau_gorx = 1, tau_res = 1)
)

nim_int_slp_samples <- runMCMC(compiled_nim_int_slp_mcmc,
                               inits = nimble_inits, nchains = 2,
                               niter = m, nburnin = bi, WAIC = TRUE)

# Print WAIC
nim_int_slp_WAIC <- nim_int_slp_samples$WAIC
print(nim_int_slp_WAIC)




bayesplot::mcmc_acf(nim_int_slp_samples$samples, pars=c('betaIncome'))

bayesplot::mcmc_dens_overlay(nim_int_slp_samples$samples)

print(MCMCvis::MCMCsummary(nim_int_slp_samples$samples))


### Forest Plot ###

library(ggplot2)

# Data preparation
# Manually inputting  from the summary output
intercepts <- data.frame(
  Group = 1:12,
  Mean = c(-0.030750479, -0.045732699, -0.059827096, -0.005951274, -0.045103576,
           0.016091896, 0.003968394, 0.025314524, 0.068392770, -0.055835231,
           -0.007908294, 0.081495252),
  SD = c(0.06233841, 0.05996935, 0.05966781, 0.06025181, 0.06024574,
         0.05987845, 0.06001767, 0.05865994, 0.06017529, 0.06283331,
         0.06075337, 0.06584236)
)

slopes <- data.frame(
  Group = 1:12,
  Mean = c(0.062307496, 0.012907488, -0.013545923, 0.026723326, -0.020827536,
           0.003419547, -0.002646676, 0.043033854, -0.048506243, -0.002535501,
           0.042886820, -0.034616441),
  SD = c(0.07346758, 0.06926529, 0.07025541, 0.07094924, 0.06915307,
         0.06938453, 0.06956820, 0.06865942, 0.07114884, 0.07420890,
         0.07087267, 0.07848699)
)

# Creating Forest plots
# Plot for intercepts
ggplot(intercepts, aes(x = Mean, y = as.factor(Group))) +
  geom_point(color = "blue") +  # Plot mean points
  geom_errorbar(aes(xmin = Mean - SD, xmax = Mean + SD), color = "blue", width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted") + # Reference line at x = 0
  labs(title = "Forest Plot for Random Intercepts",
       x = "Intercept Estimates",
       y = "Region") +
  theme_minimal()

# Plot for slopes
ggplot(slopes, aes(x = Mean, y = as.factor(Group))) +
  geom_point(color = "red") +  # Plot mean points
  geom_errorbar(aes(xmin = Mean - SD, xmax = Mean + SD), color = "red", width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted") + # Reference line at x = 0
  labs(title = "Forest Plot for Random Slopes",
       x = "Slope Estimates",
       y = "Region") +
  theme_minimal()


