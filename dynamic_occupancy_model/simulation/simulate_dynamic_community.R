library(tidyverse)

##------------------------------------------------------------------------------
# 4.3 Simulation and Analysis of dynocc model

## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_sites <- 100 # number of sites
n_years <- 7 # number of years
n_visits <- 6 # number of surveys per year

# set parameter values
psi1 <- 0.7 # prob of initial occupancy
phi <- 0.8 # persistence probability
gamma <- 0.05 # colonization probability
p0 <- -1.5 # probability of detection (logit scaled)
p_habitat_type <- 0.5 # increase in detection rate moving from one habitat type to the other (logit scaled)

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- FALSE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
  n_sites, n_years, n_visits,
  psi1, phi, gamma, 
  p0, p_habitat_type,
  create_missing_data, prob_missing
  ){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  # choose sample sizes
  n_sites <- n_sites # number of sites
  n_years <- n_years # number of years
  n_visits <- n_visits # number of surveys per year
  
  habitat_type <- rep(c(0,1), each = n_sites / 2)
  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_sites, n_years, n_visits)) # observed data
  
  # set parameter values
  psi1 <- psi1 # prob of initial occupancy
  phi <- phi # persistence probability
  gamma <- gamma # colonization probability
  #p <- p # probability of detection
  (psi_eq <- gamma / (gamma+(1-phi))) # equilibrium occupancy rate
  
  # p matrix (variable detection rate)
  logit_p_matrix <- array(NA, dim =c(n_sites, n_years))
  
  # heterogeneity in parameters
  
    for(site in 1:n_sites) { # for each interval
      for(year in 1:n_years) { # for each species
        #for(visit in 1:n_visits) { # for each visit (but sim constant rates across visits)
          
          logit_p_matrix[site, year] <- # detection is equal to 
            p0 + # an intercept
            p_habitat_type * habitat_type[site] # a spatial detection effect

        #} # for each visit
      } # for each year
    } # for each site
  
  # generate initial presence/absence states
  z[,1] <- rbinom(n=n_sites, size=1, prob=psi1) #
  sum(z[,1]) / n_sites # true occupancy proportion in year 1
  
  # generate presence/absence in subsequent years
  for(k in 2:n_years){
    
    # use z as a switch so we are estimating 
    exp_z <- z[,k-1] * phi + # survival if z=1
      (1 - z[,k-1]) * gamma # or colonization if z=0
    
    # and then assign z stochastically
    # some sites may transition if they are colonized or local extinction occurs
    # but might otherwise retain their state across years
    z[,k] <- rbinom(n = n_sites, size = 1, prob = exp_z) 
  }
  
  apply(z, 2, sum) / n_sites # true occupancy proportions in each year

  # detection / non-detection data
  for(j in 1:n_sites){
    for(k in 1:n_years){
      for(l in 1:n_visits){
        y[j,k,l] <- rbinom(n = 1, size = 1, prob = z[j,k]*ilogit(logit_p_matrix[j,k]))
      }
    }
  }

  y ; str(y)
  
  sum((y/n_visits) / sum(z)) # proportion of times detection given presence
  
  if(create_missing_data == TRUE){
    # generate missing values: create simple version of unbalanced data set
    prob_missing <- prob_missing # constant NA probability
    y2 <- y # duplicate balanced dataset
    for(i in 1:n_sites){
      for(j in 1:n_visits){
        for(t in 1:n_years){
          turnNA <- rbinom(1,1,prob_missing)
          y2[i,t,j] <- ifelse(turnNA==1, NA, y2[i,t,j])
        }
      }
    }
    y2 ; str(y2)
  } else {
    y2 <- y # duplicate dataset with no missing data
  }
  
  
  table(nvisits <- apply(y2, c(1,3), function(x) sum(!is.na(x))))
  # _ sites with no visits, _ with 1 visit, and _ with 2 visits through 7 visits
  
  # compute true expected and realized occupancy (psi and psi_fs)
  psi <- numeric(n_years) ; psi[1] <- psi1
  for(t in 2:n_years){ # compute true values of psi
    psi[t] <- psi[t-1] * phi + (1 - psi[t-1]) * gamma
  }
  psi_fs <- colSums(z) / n_sites # psi finite sample
  
  # compute observed occupancy proportion
  zobs <- apply(y2, c(1,2), function(x) max(x, na.rm = TRUE))
  zobs[zobs == "-Inf"] <- NA
  psi_obs <- apply(zobs, 2, sum, na.rm = TRUE) / apply(zobs, 2, function(x) sum(!is.na(x)))
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    psi = psi,
    psi_fs = psi_fs,
    psi_obs = psi_obs,
    V = y2, # return detection data after potentially introducing NAs,
    habitat_type = habitat_type
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

set.seed(1)
my_simulated_data <- simulate_data(  
  n_sites, n_years, n_visits,
  psi1, phi, gamma, 
  p0, p_habitat_type,
  create_missing_data, prob_missing)

V <- my_simulated_data$V
habitat_type <- my_simulated_data$habitat_type
psi <- my_simulated_data$psi
psi_fs <- my_simulated_data$psi_fs
psi_obs <- my_simulated_data$psi_obs

# Plot trajectories of psi, psi_fs, psi_eq and psi_obs
plot(1:n_years, psi, type = 'l', col = "red", lwd = 3, ylim = c(0,1), xlab = "Year",
     ylab = "Occupancy of Species speciosa", frame = F)
points(1:n_years, psi_fs, type = 'b', col = "red", pch = 16, cex = 2)
abline(h = gamma / (gamma+(1-phi)), lty = 2)
lines(1:n_years, psi_obs, col = "blue", lwd = 3)
legend('topright', c('True psi', 'True finite-sample psi', 'Equilibrium psi', 'Observed psi'),
       col = c('red', 'red', 'black', 'blue'), pch = c(NA, 16, NA, NA), lty = c(1,1,2,1),
       lwd = c(3,1,2,3), cex = 1.2, bty = 'n')


## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "habitat_type",
               "n_sites", "n_years", "n_visits") 

## Parameters monitored
params <- c("psi", "phi", "gamma", "p0", "p_habitat_type", "n_occ", "growthr", "turnover")

# MCMC settings
n_iterations <- 600
n_thin <- 1
n_burnin <- 300
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  psi1, "NA: is a vector", phi, (1-phi), gamma, p
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1 = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       p_habitat_type = runif(1, -1, 1)
  )
)

targets <- as.data.frame(cbind(params, parameter_values))


## --------------------------------------------------
### Run model

library(rstan)
stan_model <- "./dynamic_occupancy_model/models/dynocc_model0.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

print(stan_out_sim, digits = 3)
saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "psi", "phi", "gamma", "p0", "n_occ", "growthr", "turnover"
  ))




# run simpler model
## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "habitat_type",
               "n_sites", "n_years", "n_visits") 

## Parameters monitored
params <- c("psi", "phi", "gamma", "p0", "p_habitat_type", "n_occ", "growthr", "turnover")

# MCMC settings
n_iterations <- 600
n_thin <- 1
n_burnin <- 300
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  psi1, "NA: is a vector", phi, (1-phi), gamma, p
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1 = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       p_habitat_type = runif(1, -1, 1)
  )
)

targets <- as.data.frame(cbind(params, parameter_values))


## --------------------------------------------------
### Run model

library(rstan)
stan_model <- "./dynamic_occupancy_model/models/dynocc_model_simplest.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

print(stan_out_sim, digits = 3)
saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "psi", "phi", "gamma", "p0", "n_occ", "growthr", "turnover"
))