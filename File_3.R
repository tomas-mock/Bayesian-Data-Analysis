setwd("/Users/tomasmock/LSE/Bayesian Data Analysis/Reading week project")

livingexp.dat<-read.csv("LES_RW.csv",header = TRUE, stringsAsFactors = TRUE)



library(gridExtra)
library(ggplot2)
library(tidyverse)
library(arm)
library(nimble)
library(rstan)
library(rstanarm)

print(head(livingexp.dat))

livingexp.dat <- dplyr::select(livingexp.dat, c("expenditure", "income", "Gorx", "G019r", "A093r")) %>%
  mutate(G019r = factor(G019r)) %>%
  mutate(A093r = factor(A093r)) %>%
  mutate(income = ifelse(income == 0, 0.5, income)) %>%
  mutate(expenditure = log(expenditure)) %>%
  mutate(income = log(income))

ggplot(livingexp.dat, aes(x=income,y=expenditure,color=G019r)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + facet_wrap(~Gorx)

livingexp.dat %>% filter(income>2.5) %>% ggplot(aes(x=income,y=expenditure,color=G019r)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + facet_wrap(~Gorx)

no_low_livingexp.dat <- livingexp.dat %>% filter(income>2.5)



# centre the expenditure and income variables
livingexp.dat.tr <- no_low_livingexp.dat %>%
  mutate(expenditure = expenditure - mean(expenditure)) %>%
  mutate(income = income - mean(income))

# Read and prepare your data
livingexp.dat.tr <- livingexp.dat.tr %>%
  mutate(A093r = relevel(A093r, "1"), G019r = relevel(G019r, "1"))

# Model 2 with STAN #

# Create model matrix for fixed effects
X <- model.matrix(expenditure ~ income + A093r + G019r, data = livingexp.dat.tr)

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
  list(beta = rep(0, D), nu_intercept = rep(0, nGorx), sigma_res = 1, sigma_nu_intercept = 1),
  list(beta = rep(1, D), nu_intercept = rep(1, nGorx), sigma_res = 10, sigma_nu_intercept = 10)
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
  real<lower=0> sigma_nu_intercept;
}

model {
  beta ~ normal(0, 2.5);
  nu_intercept ~ normal(0, sigma_nu_intercept);
  sigma_res ~ normal(0, 5);
  sigma_nu_intercept ~ normal(0, 5);
  for (i in 1:N) {
    expenditure[i] ~ normal(dot_product(beta, X[i]) + nu_intercept[Gorx[i]], sigma_res);
  }
}

generated quantities {
  vector[N] expenditure_pred;
  vector[N] log_lik;
  for (i in 1:N) {
    expenditure_pred[i] = normal_rng(dot_product(beta, X[i]) + nu_intercept[Gorx[i]], sigma_res);
    log_lik[i] = normal_lpdf(expenditure[i] | dot_product(beta, X[i]) + nu_intercept[Gorx[i]], sigma_res);
  }
}
")

# Compile and sample from the model
hierarchical_stan_model <- stan_model("stan_models/hierarchical_living.stan", model_name = "hierarchical_living")
hierarchical_living_stan_fit <- sampling(hierarchical_stan_model, data = hierarchical_stan_data, warmup = 4000, iter = 10000, chains = 2, thin = 2, init = hierarchical_stan_inits)


### summary ###
fit_summary <- summary(hierarchical_living_stan_fit)

print(fit_summary)






# Extract predictions
expenditure_pred_extracted2 <- extract(hierarchical_living_stan_fit)$expenditure_pred

# visualize the entire distribution of predictions
ppc_dens_overlay(y = livingexp.dat.tr$expenditure, yrep = expenditure_pred_extracted2)






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





# More plots #

bayesplot::mcmc_pairs(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]', 'beta[5]', 'beta[6]', 'beta[7]'))


bayesplot::mcmc_acf(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]', 'beta[5]', 'beta[6]', 'beta[7]'))


bayesplot::mcmc_trace(hierarchical_living_stan_fit, pars=c('beta[1]','beta[2]','beta[3]','beta[4]', 'beta[5]', 'beta[6]', 'beta[7]'))



# MODEL 2 with NIMBLE # 
# Convert Gorx 
livingexp.dat.tr$Gorx <- as.factor(livingexp.dat.tr$Gorx)


summary(livingexp.dat.tr)

# Model setup
nGroups <- length(unique(livingexp.dat.tr$Gorx)) # Number of Gorx groups
N <- nrow(livingexp.dat.tr)                      # Number of observations
J <- length(unique(livingexp.dat.tr$A093r))      # Number of levels in A093r
K <- length(unique(livingexp.dat.tr$G019r))      # Number of levels in G019r

# Create a NIMBLE model code block
nim_int_slp <- nimbleCode({
  beta0 ~ dnorm(0, 0.25) 
  betaIncome ~ dnorm(0, 0.25)
  for (j in 1:J) {
    betaA093r[j] ~ dnorm(0, 0.25)
  }
  for (k in 1:K) {
    betaG019r[k] ~ dnorm(0, 0.25)
  }
  for (g in 1:nGroups) {
    interceptGorx[g] ~ dnorm(0, tau_gorx)
  }
  tau_gorx ~ dgamma(2.5, 0.5)
  sigma_gorx <- sqrt(1 / tau_gorx)
  tau_res ~ dgamma(2.5, 0.5)
  sigma_res <- sqrt(1 / tau_res)
  for (i in 1:N) {
    mu[i] <- beta0 + betaIncome * income[i] +
      inprod(betaA093r[1:J], A093r[i,1:J]) +
      inprod(betaG019r[1:K], G019r[i,1:K]) +
      interceptGorx[GorxID[i]]
    expenditure[i] ~ dnorm(mu[i], tau_res)
  }
})

# Prepare data and initial values for NIMBLE
A093r_matrix <- model.matrix(~ A093r - 1, data = livingexp.dat.tr)
G019r_matrix <- model.matrix(~ G019r - 1, data = livingexp.dat.tr)
Gorx_indices <- as.numeric(livingexp.dat.tr$Gorx)

nimble_data <- list(expenditure = livingexp.dat.tr$expenditure,
                    income = livingexp.dat.tr$income,
                    A093r = A093r_matrix,
                    G019r = G019r_matrix,
                    GorxID = Gorx_indices)

nimble_constants <- list(N = N, J = J, K = K, nGroups = nGroups)
nimble_inits <- list(beta0 = 0, betaIncome = 0, betaA093r = rep(0, J),
                     betaG019r = rep(0, K), interceptGorx = rep(0, nGroups),
                     tau_gorx = 1, tau_res = 1)

# Build and compile the model with WAIC enabled
nim_int_slp_model <- nimbleModel(code = nim_int_slp, 
                                 constants = nimble_constants, 
                                 data = nimble_data, 
                                 inits = nimble_inits)
compiled_nim_int_slp_model <- compileNimble(nim_int_slp_model)
nim_int_slp_mcmc <- buildMCMC(nim_int_slp_model, enableWAIC = TRUE)
compiled_nim_int_slp_mcmc <- compileNimble(nim_int_slp_mcmc, project = nim_int_slp_model)

# Run the MCMC with WAIC calculation
nimble_inits_multiple_chains <- list(
  list(beta0 = 0, betaIncome = 0.1, betaA093r = rep(0.1, J), betaG019r = rep(0.1, K), 
       interceptGorx = rep(0.1, nGroups), tau_gorx = 1, tau_res = 1),
  list(beta0 = 0, betaIncome = -0.1, betaA093r = rep(-0.1, J), betaG019r = rep(-0.1, K), 
       interceptGorx = rep(-0.1, nGroups), tau_gorx = 1, tau_res = 1)
)

nim_int_slp_samples <- runMCMC(compiled_nim_int_slp_mcmc,
                               inits = nimble_inits_multiple_chains, nchains = 2,
                               niter = m, nburnin = bi, WAIC = TRUE)

# Print WAIC
nim_int_slp_WAIC <- nim_int_slp_samples$WAIC
print(nim_int_slp_WAIC)



bayesplot::mcmc_acf(nim_int_slp_samples$samples, pars=c('betaIncome'))

bayesplot::mcmc_dens_overlay(nim_int_slp_samples$samples)

print(MCMCvis::MCMCsummary(nim_int_slp_samples$samples))

# Data preparation
intercepts <- data.frame(
  Group = 1:12,
  Mean = c(-0.064093801, -0.055564135, -0.064067588, -0.010890015, -0.052918336,
           0.018416300, -0.006727243, 0.033895242, 0.063114759, -0.059601361,
           -0.019806858, 0.080673001),
  SD = c(0.07932518, 0.07745362, 0.07736156, 0.07824788, 0.07769467,
         0.07723161, 0.07788650, 0.07721846, 0.07837932, 0.07928236,
         0.07803521, 0.08117488)
)

# Creating Forest Plots plots
# Plot for intercepts
ggplot(intercepts, aes(x = Mean, y = as.factor(Group))) +
  geom_point(color = "blue") +  # Plot mean points
  geom_errorbar(aes(xmin = Mean - SD, xmax = Mean + SD), color = "blue", width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted") + # Reference line at x = 0
  labs(title = "Forest Plot for Random Intercepts",
       x = "Intercept Estimates",
       y = "Group") +
  theme_minimal()



