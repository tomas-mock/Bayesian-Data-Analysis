setwd("/Users/tomasmock/LSE/Bayesian Data Analysis/Reading week project")

livingexp.dat<-read.csv("LES_RW.csv",header = TRUE, stringsAsFactors = TRUE)

print(head(livingexp.dat))

library(gridExtra)
library(ggplot2)
library(tidyverse)
library(arm)
library(nimble)
library(rstan)
library(rstanarm)
library(MCMCvis)



livingexp.dat <- dplyr::select(livingexp.dat, c("expenditure","income","Gorx","A121r", "A093r", "G019r"))%>% mutate(A121r=factor(A121r))%>%mutate(A093r=factor(A093r))%>%mutate(G019r=factor(G019r))%>%
  mutate(income=ifelse(income==0,0.5,income))%>%
  mutate(expenditure=log(expenditure)) %>%
  mutate(income=log(income))

ggplot(livingexp.dat, aes(x=income,y=expenditure,color=A121r)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + facet_wrap(~Gorx)

livingexp.dat %>% filter(income>2.5) %>% ggplot(aes(x=income,y=expenditure,color=A121r)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + facet_wrap(~Gorx)

no_low_livingexp.dat <- livingexp.dat %>% filter(income>2.5)


# Exploratory data analysis

print(head(no_low_livingexp.dat))

summary(no_low_livingexp.dat)

str(no_low_livingexp.dat)

# Histogram of log Expenditure
ggplot(data = no_low_livingexp.dat, aes(x = expenditure)) + 
  geom_histogram(binwidth = 0.2, fill = "#0000ff97", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of log Expenditure", x = "log Expenditure", y = "Frequency")

# Histogram of log Income
ggplot(data = no_low_livingexp.dat, aes(x = income)) + 
  geom_histogram(binwidth = 0.2, fill = "#0000ffb4", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of log Income", x = "log Income", y = "Frequency")

# Bar chart of 'A121r'
ggplot(data = no_low_livingexp.dat, aes(x = A121r)) + 
  geom_bar() +
  theme_minimal() +
  labs(title = "Bar Chart of A121r", x = "A121r", y = "Count")

# Bar chart of 'A093r'
ggplot(data = no_low_livingexp.dat, aes(x = A093r)) + 
  geom_bar() +
  theme_minimal() +
  labs(title = "Bar Chart of A093r", x = "A093r", y = "Count")

# Bar chart of 'G019r'
ggplot(data = no_low_livingexp.dat, aes(x = G019r)) + 
  geom_bar() +
  theme_minimal() +
  labs(title = "Bar Chart of G019r", x = "G019r", y = "Count")


# Bar chart of 'Gorx'
ggplot(data = no_low_livingexp.dat, aes(x = Gorx)) + 
  geom_bar() +
  theme_minimal() +
  labs(title = "Bar Chart of Gorx", x = "Gorx", y = "Count")

# Scatterplot of 'log Expenditure' vs 'log Income'
ggplot(data = no_low_livingexp.dat, aes(x = income, y = expenditure)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "Scatterplot of log Expenditure vs log Income", x = "log Income", y = "log Expenditure")

# Comparative boxplot of 'log Expenditure' for each level of 'A121r'
ggplot(data = no_low_livingexp.dat, aes(x = A121r, y = expenditure)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparative Boxplot of log Expenditure by A121r", x = "A121r", y = "log Expenditure")


# Comparative boxplot of 'log Expenditure' for each level of 'A093r'
ggplot(data = no_low_livingexp.dat, aes(x = A093r, y = expenditure)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparative Boxplot of log Expenditure by A093r", x = "A093r", y = "log Expenditure")

# Comparative boxplot of 'log Expenditure' for each level of 'G019r'
ggplot(data = no_low_livingexp.dat, aes(x = G019r, y = expenditure)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparative Boxplot of log Expenditure by G019r", x = "G019r", y = "log Expenditure")



# MODEL 1-F  FREQUENTIST LINEAR REGRESSION
# Fit the linear model
model <- lm(expenditure ~ income + A121r, data = no_low_livingexp.dat)

# Summarize the model to see results
summary(model)

# 1. Residuals vs Fitted
plot(model, which = 1)

# 2. Normal Q-Q
plot(model, which = 2)

# 3. Scale-Location (Spread-Location)
plot(model, which = 3)

# 4. Residuals vs Leverage
plot(model, which = 4)


# MODEL 2-F  FREQUENTIST LINEAR REGRESSION
# Fit the linear model
model <- lm(expenditure ~ income + A093r + G019r, data = no_low_livingexp.dat)

# Summarize the model to see results
summary(model)

# 1. Residuals vs Fitted
plot(model, which = 1)

# 2. Normal Q-Q
plot(model, which = 2)

# 3. Scale-Location (Spread-Location)
plot(model, which = 3)

# 4. Residuals vs Leverage
plot(model, which = 4)



# centre the age and income variables
livingexp.dat.tr <- no_low_livingexp.dat %>%
  mutate(expenditure = expenditure - mean(expenditure)) %>%
  mutate(income = income - mean(income))


### MODEL 1-N (BAYESIAN REGRESSION WITH NIMBLE) ###

# Create dummies for categorical variables
X <- model.matrix(~ income + A121r, data = livingexp.dat.tr)
print(head(X))

living_nim_code <- nimbleCode({
  # Priors
  # Regression coefficients
  for (d in 1:D){
    beta[d] ~ dnorm(0, 0.0001)
  }
  # Residual variance
  tau ~ dgamma(2.5, 0.5)
  sigma <- sqrt(1/tau)
  # Likelihood
  for (i in 1:N){
    mu[i] <- beta[1] * X[i, 1] + beta[2] * X[i, 2]
    expenditure[i] ~ dnorm(mu[i], tau)
  }
})

nimble_constants <- list(N = nrow(X), X = X, D = ncol(X))
nimble_data <- list(expenditure = livingexp.dat.tr$expenditure)
nimble_inits <- list(beta = rep(0, ncol(X)), tau = 0.1)

# Create and compile the NIMBLE model
living_nim_model <- nimbleModel(code = living_nim_code, name = "living_nimble_project",
                                constants = nimble_constants, data = nimble_data, inits = nimble_inits) 

# Recompile the model with resetFunctions set to TRUE
compiled_liv_nim_model <- compileNimble(living_nim_model, resetFunctions = TRUE)

# Build the MCMC object with WAIC enabled
liv_nim_mcmc <- buildMCMC(compiled_liv_nim_model, enableWAIC = TRUE)

# Recompile the MCMC with the updated compiled model
compiled_liv_nim_mcmc <- compileNimble(liv_nim_mcmc, project = compiled_liv_nim_model)

# Set up multiple initial conditions for MCMC
nimble_inits = list(list(beta = rep(0, ncol(X)), tau = 0.1), list(beta = rep(1, ncol(X)), tau = 1)) 

# MCMC settings
m <- 10000
bi <- 4000

# Run MCMC with WAIC calculation
liv_nim_samples <- runMCMC(compiled_liv_nim_mcmc,
                           inits = nimble_inits, nchains = 2,
                           niter = m, nburnin = bi, WAIC = TRUE)

# Retrieve WAIC and print it
liv_nim_WAIC <- liv_nim_samples$WAIC
print(liv_nim_WAIC)


library(MCMCvis)

print(MCMCvis::MCMCsummary(liv_nim_samples$samples))

beta_with <- MCMCvis::MCMCsummary(liv_nim_samples$samples)[1:2,1]

bayesplot::mcmc_trace(liv_nim_samples$samples,
                      facet_args = list(ncol = 1, strip.position = "left"))
#density plot
bayesplot::mcmc_dens_overlay(liv_nim_samples$samples) 
#paired plots
bayesplot::mcmc_pairs(liv_nim_samples$samples) 
#autocorrelation functions
bayesplot::mcmc_acf(liv_nim_samples$samples)
#gelman diagnostic
coda::gelman.plot(liv_nim_samples$samples)





### MODEL 1-S (BAYESIAN REGRESSION WITH STAN) ###

livingexp.dat.tr <- livingexp.dat.tr %>%
  mutate(A121r = relevel(A121r, "1"))


X <- model.matrix(expenditure ~ income + A121r, data = livingexp.dat.tr)

living_stan_data <- list(expenditure=livingexp.dat.tr$expenditure, X=X, N=nrow(X), D=ncol(X))
living_stan_inits <- list(list(beta=rep(0, ncol(X)), sigma=1),list(beta=rep(10,ncol(X), sigma=10)))

print(head(X))

# Creating the directory for the Stan models
if (!dir.exists("stan_models")) {
  dir.create("stan_models")
}


write("
// Stan model for linear regression
data {
  int<lower=1> N;   // Number of data items
  int<lower=1> D;   // Number of predictors
  matrix[N, D] X;   // Predictor matrix
  vector[N] expenditure; // Outcome vector
}

parameters {
  vector[D] beta; // Coefficients for predictors
  real<lower=0> sigma; // SD of residuals
}

model {
  // Priors
  beta ~ normal(0, 2.5);
  sigma ~ normal(0, 5); // 

  // Likelihood
  expenditure ~ normal(X * beta, sigma);
}
", file = "stan_models/living.stan")


# Compiling the Stan model
living_stan_model <- stan_model("stan_models/living.stan", model_name = "living_stan")

living_stan_fit <- sampling(living_stan_model, data = living_stan_data,
                            warmup = 4000,
                            iter = 10000,
                            chains = 2,
                            thin = 2,
                            init = living_stan_inits)

bayesplot::mcmc_pairs(living_stan_fit) 
bayesplot::mcmc_acf(living_stan_fit) 
bayesplot::mcmc_trace(living_stan_fit)

library(coda)
living_stan_hmc <- As.mcmc.list(living_stan_fit)
gelman.diag(living_stan_hmc) 
gelman.plot(living_stan_hmc) 

output <- as.data.frame(summary(living_stan_fit, pars=c("sigma",paste("beta[",1:4,"]",sep="")), probs=c(0.025, 0.975))$summary)



rownames(output)<-c("sigma",colnames(X)) 


print(output)




### MODEL 2-N (BAYESIAN REGRESSION WITH NIMBLE) ###

# Update the model matrix to include new predictors and exclude the old one
X <- model.matrix(~ income + A093r + G019r, data = livingexp.dat.tr)

# Display the first few rows of the model matrix
print(head(X))

# Define the model structure using nimbleCode
living_nim_code <- nimbleCode({
  # Priors and model definition
  for (d in 1:D) {
    beta[d] ~ dnorm(0, 0.0001)  # priors for regression coefficients
  }
  tau ~ dgamma(2.5, 0.5)  # prior for residual variance
  sigma <- sqrt(1/tau)  # standard deviation
  # Likelihood calculation
  for (i in 1:N) {
    mu[i] <- beta[1] * X[i, 1] + beta[2] * X[i, 2] + beta[3] * X[i, 3]
    expenditure[i] ~ dnorm(mu[i], tau)
  }
})

# Constants, data, and initial values for the NIMBLE model
nimble_constants <- list(N = nrow(X), X = X, D = ncol(X))
nimble_data <- list(expenditure = livingexp.dat.tr$expenditure)
nimble_inits <- list(beta = rep(0, ncol(X)), tau = 0.1)

# Create and compile the NIMBLE model
living_nim_model <- nimbleModel(code = living_nim_code, name = "living_nimble_project", 
                                constants = nimble_constants, data = nimble_data, inits = nimble_inits)

# Recompile the model with resetFunctions set to TRUE
compiled_liv_nim_model <- compileNimble(living_nim_model, resetFunctions = TRUE)

# Build the MCMC object with WAIC enabled
liv_nim_mcmc <- buildMCMC(compiled_liv_nim_model, enableWAIC = TRUE)

# Recompile the MCMC with the updated compiled model
compiled_liv_nim_mcmc <- compileNimble(liv_nim_mcmc, project = compiled_liv_nim_model)

# Set up multiple initial conditions for MCMC
nimble_inits <- list(list(beta = rep(0, ncol(X)), tau = 0.1), list(beta = rep(1, ncol(X)), tau = 1))

# MCMC settings
m <- 10000
bi <- 4000

# Run MCMC with WAIC calculation
liv_nim_samples <- runMCMC(compiled_liv_nim_mcmc,
                           inits = nimble_inits, nchains = 2,
                           niter = m, nburnin = bi, WAIC = TRUE)

# Retrieve WAIC and print it
liv_nim_WAIC <- liv_nim_samples$WAIC
print(liv_nim_WAIC)



print(MCMCvis::MCMCsummary(liv_nim_samples$samples))

beta_with <- MCMCvis::MCMCsummary(liv_nim_samples$samples)[1:2,1]

bayesplot::mcmc_trace(liv_nim_samples$samples,
                      facet_args = list(ncol = 1, strip.position = "left"))
#density plot
bayesplot::mcmc_dens_overlay(liv_nim_samples$samples) 
#paired plots
bayesplot::mcmc_pairs(liv_nim_samples$samples) 
#autocorrelation functions
bayesplot::mcmc_acf(liv_nim_samples$samples)
#gelman diagnostic
coda::gelman.plot(liv_nim_samples$samples)





### MODEL 2-S (Bayesian Regression with STAN) ###


# A093r and G019r as factors
livingexp.dat.tr <- livingexp.dat.tr %>%
  mutate(A093r = factor(A093r), G019r = factor(G019r))

# Update model matrix to include new predictors as factors
X <- model.matrix(~ income + A093r + G019r, data = livingexp.dat.tr)


# Prepare data for Stan
living_stan_data <- list(expenditure = livingexp.dat.tr$expenditure, X = X, N = nrow(X), D = ncol(X))

# Set initial values for the parameters
living_stan_inits <- list(list(beta = rep(0, ncol(X)), sigma = 1), list(beta = rep(10, ncol(X)), sigma = 10))

# Print the first few rows of the updated model matrix
print(head(X))

# Check and create a directory for the Stan models if it doesn't exist
if (!dir.exists("stan_models")) {
  dir.create("stan_models")
}

#  Stan model for linear regression
write("
// Stan model for linear regression
data {
  int<lower=1> N;   // Number of data items
  int<lower=1> D;   // Number of predictors
  matrix[N, D] X;   // Predictor matrix
  vector[N] expenditure; // Outcome vector
}

parameters {
  vector[D] beta; // Coefficients for predictors
  real<lower=0> sigma; // SD of residuals
}

model {
  // Priors
  beta ~ normal(0, 2.5);
  sigma ~ normal(0, 5); // Careful with this prior for sigma since it allows negative values

  // Likelihood
  expenditure ~ normal(X * beta, sigma);
}
", file = "stan_models/living.stan")

# Compile the Stan model
living_stan_model <- stan_model("stan_models/living.stan", model_name = "living_stan")

# Execute sampling
living_stan_fit <- sampling(living_stan_model, data = living_stan_data,
                            warmup = 4000,
                            iter = 10000,
                            chains = 2,
                            thin = 2,
                            init = living_stan_inits)

bayesplot::mcmc_pairs(living_stan_fit) 
bayesplot::mcmc_acf(living_stan_fit) 
bayesplot::mcmc_trace(living_stan_fit)

library(coda)
living_stan_hmc <- As.mcmc.list(living_stan_fit)
gelman.diag(living_stan_hmc) 
gelman.plot(living_stan_hmc) 

output <- as.data.frame(summary(living_stan_fit, pars=c("sigma",paste("beta[",1:7,"]",sep="")), probs=c(0.025, 0.975))$summary)



rownames(output)<-c("sigma",colnames(X)) 


print(output)




