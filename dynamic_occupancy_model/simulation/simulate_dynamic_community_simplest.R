
##------------------------------------------------------------------------------
# 4.3 Simulation and Analysis of the simplest dynocc model

## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_sites <- 150 # number of sites
n_years <- 7 # number of years
n_visits <- 6 # number of surveys per year

# set parameter values
psi1 <- 0.7 # prob of initial occupancy
phi0 <- 0.75 # persistence probability
phi_habitat <- 0.15 # persistence probability
gamma <- 0.05 # colonization probability

p0 <- -1.75 # probability of detection (logit scaled)
p_habitat <- 0.5 # increase in detection rate moving from one habitat type to the other (logit scaled)
p_date <- 0
p_date_sq <- -0.5

mean_survey_date <- 180
sigma_survey_date <- 40

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- FALSE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
    n_sites, n_years, n_visits,
    psi1, phi0, phi_habitat, gamma, 
    p0, p_habitat_type,
    p_date, p_date_sq,
    mean_survey_date, sigma_survey_date,
    create_missing_data, prob_missing
){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  # choose sample sizes
  n_sites <- n_sites # number of sites
  n_years <- n_years # number of years
  n_visits <- n_visits # number of surveys per year

  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_sites, n_years, n_visits)) # observed data
  
  # set parameter values
  psi1 <- psi1 # prob of initial occupancy
  phi0 <- phi0 # persistence probability
  phi_habitat <- phi_habitat # effect of habitat type on persistence
  gamma <- gamma # colonization probability
  # p <- p # probability of detection
  (psi_eq <- gamma / (gamma+(1-phi))) # equilibrium occupancy rate
  
  ## --------------------------------------------------
  ## Create covariate data
  
  ## day of year
  # should have a value for each site*year*visit
  date <- array(NA, dim =c(n_sites, n_years, n_visits))
  mean_survey_date = mean_survey_date
  sigma_survey_date = sigma_survey_date
  
  for(site in 1:n_sites) { # for each site
    for(n_years in 1:n_years) { # for each n_years
      # create a vector of visit dates centered on the middle of the early summer
      date[site, n_years, 1:n_visits] <- sort(as.integer(rnorm(
        n_visits, mean = mean_survey_date, sd = sigma_survey_date))) 
    }
  }
  
  ## let's scale the calendar day of year by the mean date (z-score scaled)
  date_scaled <- center_scale(date) 
  
  ## habitat type
  habitat_type <- rep(c(0,1), each = n_sites / 2)
  
  # generate p with heterogeneity
  logit_p <- array(NA, dim = c(n_sites, n_years, n_visits)) 
  
  for(j in 1:n_sites){
    for(k in 1:n_years){
      for(l in 1:n_visits){
        
        logit_p[j,k,l] = p0 +
          p_habitat * habitat_type[j] +
          p_date*date_scaled[j, k, l] + # a spatiotemporally specific intercept
          p_date_sq*(date_scaled[j, k, l])^2 # a spatiotemporally specific intercept
        
      }
    }
  }
  
  # generate initial presence/absence states
  z[,1] <- rbinom(n=n_sites, size=1, prob=psi1) #
  sum(z[,1]) / n_sites # true occupancy proportion in year 1
  
  # generate p with heterogeneity
  logit_phi <- array(NA, dim = c(n_sites, n_years)) 
  
  for(j in 1:n_sites){
    for(k in 1:n_years){

        logit_pphi[j,k] = phi0 +
          phi_habitat * habitat_type[j] 
        
    }
  }
  
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
        y[j,k,l] <- rbinom(n = 1, size = 1, prob = z[j,k]*ilogit(logit_p[j,k,l]))
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
    habitat_type = habitat_type,
    date_scaled = date_scaled
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

set.seed(1)
my_simulated_data <- simulate_data(  
  n_sites, n_years, n_visits,
  psi1, phi0, phi_habitat, gamma, 
  p0, p_habitat_type,
  p_date, p_date_sq,
  mean_survey_date, sigma_survey_date,
  create_missing_data, prob_missing)

V <- my_simulated_data$V
habitat_type <- my_simulated_data$habitat_type
psi <- my_simulated_data$psi
psi_fs <- my_simulated_data$psi_fs
psi_obs <- my_simulated_data$psi_obs
date_scaled <- my_simulated_data$date_scaled

# Plot trajectories of psi, psi_fs, psi_eq and psi_obs
plot(1:n_years, psi, type = 'l', col = "red", lwd = 3, ylim = c(0,1), xlab = "Year",
     ylab = "Occupancy of Species speciosa", frame = F)
points(1:n_years, psi_fs, type = 'b', col = "red", pch = 16, cex = 2)
abline(h = gamma / (gamma+(1-phi)), lty = 2)
lines(1:n_years, psi_obs, col = "blue", lwd = 3)
legend('topright', c('True psi', 'True finite-sample psi', 'Equilibrium psi', 'Observed psi'),
       col = c('red', 'red', 'black', 'blue'), pch = c(NA, 16, NA, NA), lty = c(1,1,2,1),
       lwd = c(3,1,2,3), cex = 1.2, bty = 'n')


# run simpler model
## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", 
               "n_sites", "n_years", "n_visits",
               "habitat_type", "date_scaled") 

## Parameters monitored
params <- c("phi", "gamma", "psi1",
            "p0", "p_habitat", "p_date", "p_date_sq",
            "T_rep", "T_obs", "P_site")

# MCMC settings
n_iterations <- 1000
n_thin <- 1
n_burnin <- 500
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  phi, gamma, psi1, p0, p_habitat, p_date, p_date_sq, NA, NA, NA
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1 = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       p_habitat = runif(1, -1, 1),
       p_date = runif(1, -1, 1),
       p_date_sq = runif(1, -1, 1)
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

print(stan_out_sim, digits = 3, 
      pars = c("phi", "gamma", "psi1",
               "p0", "p_habitat", "p_date", "p_date_sq"))

saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "psi1", "phi", "gamma", "p0", "p_habitat"
))


### PPC's

print(stan_out_sim, digits = 3, pars = c("P_site"))

fit_summary <- rstan::summary(stan_out_sim)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[308:457,1]))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)

# Evaluation of fit 
# site 1
plot(list_of_draws[,8], list_of_draws[,158], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

# site 2
plot(list_of_draws[,9], list_of_draws[,159], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

