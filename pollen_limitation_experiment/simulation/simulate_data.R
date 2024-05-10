## Simulate data for mixed-effects, logistic regression
#  started Dec. 14, 2022, J Ulrich

# Simulate binary data where: 
# the outcome is: 
# 0 - a flower is not pollen limited, OR
# 1 - a flower is pollen limited
# A random intercept for planter pot, nested in site will account for heterogeneity
# in nested sample groups.

###-----------------------------------------------------------------------------
## Global simulation options

n_sites = 20 # number of sites
n_pots_per_site = 10 # pots per site
n_pots = n_sites*n_pots_per_site
n_obs_per_pot = 10 # pairs per pot
N = n_pots_per_site*n_sites*n_obs_per_pot # number of paired flowers
alpha0 = -0.5 # global intercept for outcome 
sigma_alpha_pot = 1 # variation across sites
sigma_alpha_site = 0.5 # variation across sites
beta = 1 # effect of site type on outcome

###-----------------------------------------------------------------------------
## Simulation function

simulate_data <- function(
  
  N = N,
  n_obs_per_pot = n_obs_per_pot,
  n_pots_per_site = n_pots_per_site,
  n_pots = n_pots,
  n_sites = n_sites,
  alpha0 = alpha0,
  sigma_alpha_pot = sigma_alpha_pot,
  sigma_alpha_site = sigma_alpha_site,
  beta = beta
  
){
  
  # Inverse logit functio
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  # Generate covariate data (site type)
  x = rep(c(0, 1), each = 0.5*N)
  
  ## site-specific random intercepts
  sites = rep(1:n_sites, each = n_pots_per_site*n_obs_per_pot)
  
  # site baseline success is drawn from a normal distribution with mean 0 and 
  # site specific variation defined by sigma_alpha_site
  site_intercepts <- rep(rnorm(n=n_sites, mean=0, sd=sigma_alpha_site),
                         each=n_pots_per_site*n_obs_per_pot)
  
  ## pot-specific random intercepts
  pots = rep(1:n_pots, each = n_obs_per_pot)
  
  # site baseline success is drawn from a normal distribution with mean 0 and 
  # site specific variation defined by sigma_alpha_site
  pot_intercepts <- rep(rnorm(n=n_pots, mean=0, sd=sigma_alpha_pot),
                        each=n_obs_per_pot)

  # View(as.data.frame(cbind(sites, pots, site_intercepts, pot_intercepts)))
  
  nested_intercept <- vector(length=N)
  
  for(i in 1:N){
      
      nested_intercept[i] <- site_intercepts[i] + pot_intercepts[i]
    
  }
  
  siteLookup <- rep(1:n_sites, each=n_pots_per_site)
  
  # Generate probability that outcome is 1
  p = vector(length=N)
  
  for(i in 1:N){
    
    p[i] = 
      alpha0 + # a global intercept
      nested_intercept[i] + # a site-specific intercept adjustment
      beta * x[i] # plus an effect of site type
    
  }
  
  
  # Generate outcome data
  y = vector(length = N)
  
  for(i in 1:N){
    
    y[i] = rbinom(1, 1, prob=inv_logit(p[i])) 
    
  }
  
  ###-----------------------------------------------------------------------------
  ## Return stuff
  
  return(list(
    
    N = N, # number of pairs
    n_pots_per_site = n_pots_per_site, # number of pots
    n_pots = n_pots,
    pots = pots, # vector of pot names
    n_sites = n_sites, # number of sites
    sites = sites, # vector of site names
    siteLookup = siteLookup,
    x = x, # site type covariate data
    y = y # outcome data
    
  ))
  
} # end simulate_data() function


## --------------------------------------------------
### Simulate data
set.seed(3)
my_simulated_data <- simulate_data(N,
                                   n_obs_per_pot,
                                   n_pots_per_site,
                                   n_pots,
                                   n_sites,
                                   alpha0,
                                   sigma_alpha_pot,
                                   sigma_alpha_site,
                                   beta)


## --------------------------------------------------
### Prepare data for model

# data to feed to the model
N <- my_simulated_data$N # number of pairs
n_pots <- my_simulated_data$n_pots # number of sites
pot_ID <- as.numeric(my_simulated_data$pots) # number of sites
n_sites <- my_simulated_data$n_sites # number of sites
site_ID <- as.numeric(my_simulated_data$sites) # vector of sites
siteLookup <- as.numeric(my_simulated_data$siteLookup) 
x <- my_simulated_data$x # site type covariate
y <- my_simulated_data$y # outcome

stan_data <- c("N", 
               "n_pots", "pot_ID",
               "n_sites", "site_ID",
               "siteLookup",
               "x", "y")

# Parameters monitored
params <- c("alpha0",
            "beta",
            "sigma_alpha_site",
            "sigma_alpha_pot"
)

parameter_value <- c(alpha0,
                     beta,
                     sigma_alpha_site,
                     sigma_alpha_pot
)

# MCMC settings
n_iterations <- 2000
n_thin <- 2
n_burnin <- 0.5*n_iterations
n_chains <- 3
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(alpha0 = runif(1, -1, 1),
       sigma_alpha_site = runif(1, 0, 1),
       sigma_alpha_pot = runif(1, 0, 1),
       beta = runif(1, -1, 1)
       
  )
)

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./models/logistic_model.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores,
                     control = list(adapt_delta = 0.9))

print(stan_out_sim, digits = 3)
View(targets)

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "alpha0",
  "beta",
  "sigma_alpha_site",
  "sigma_alpha_pot"
))

# pairs plot
pairs(stan_out_sim, pars = c(
  "alpha0",
  "beta",
  "sigma_alpha_site",
  "sigma_alpha_pot"
))

