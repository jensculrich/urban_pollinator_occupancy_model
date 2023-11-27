
## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_species <- 60 # number of species
n_sites <- 60 # number of sites
n_years <- 7 # number of years
n_visits <- 6 # number of surveys per year

# set parameter values
psi1 <- 0.7 # prob of initial occupancy
phi0 <- 1 # persistence probability
sigma_phi_species <- 1 # species-specific variation
phi_habitat <- 1.5 # persistence probability
gamma <- 0.075 # colonization probability

p0 <- 0.5 # probability of detection (logit scaled)
sigma_p_species <- 1 # species-specific variation
p_habitat <- -0.5 # increase in detection rate moving from one habitat type to the other (logit scaled)
mu_p_species_date <- 0
sigma_p_species_date <- 0.5
mu_p_species_date_sq <- -0.5  
sigma_p_species_date_sq <- 0.5

mean_survey_date <- 180
sigma_survey_date <- 40

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- FALSE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
    n_species, n_sites, n_years, n_visits,
    psi1, phi0, sigma_phi_species, phi_habitat, gamma, 
    p0, sigma_p_species, p_habitat_type,
    mu_p_species_date, sigma_p_species_date, mu_p_species_date_sq, sigma_p_species_date_sq,
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
  n_species <- n_species # number of species 
  n_sites <- n_sites # number of sites
  n_years <- n_years # number of years
  n_visits <- n_visits # number of surveys per year
  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_species, n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_species, n_sites, n_years, n_visits)) # observed data
  
  # set parameter values
  psi1 <- psi1 # prob of initial occupancy
  phi0 <- phi0 # persistence probability
  sigma_phi_species <- sigma_phi_species # species-specific variation
  phi_habitat <- phi_habitat # effect of habitat type on persistence
  gamma <- gamma # colonization probability
  # p <- p # probability of detection
  # equilibrium occupancy rate (for the average species)
  (psi_eq_habitat0 <- gamma / (gamma+(1-ilogit(phi0 + phi_habitat*0)))) # equilibrium occupancy rate 
  (psi_eq_habitat1 <- gamma / (gamma+(1-ilogit(phi0 + phi_habitat*1)))) # equilibrium occupancy rate
  
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
  

  ## --------------------------------------------------
  ## Create random effects
  
  ## species-specific random intercepts
  phi_species <- rnorm(n=n_species, mean=phi0, sd=sigma_phi_species)
  
  p_species <- rnorm(n=n_species, mean=p0, sd=sigma_p_species)
  
  p_species_date <- rnorm(n=n_species, mean=mu_p_species_date, sd=sigma_p_species_date)
  
  p_species_date_sq <- rnorm(n=n_species, mean=mu_p_species_date_sq, sd=sigma_p_species_date_sq)
  
  ## --------------------------------------------------
  ## Create expected values
  
  # generate p with heterogeneity
  logit_p <- array(NA, dim = c(n_species, n_sites, n_years, n_visits)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_visits){
          
          logit_p[i,j,k,l] = 
            p_species[i] +
            p_habitat * habitat_type[j] +
            p_species_date[i]*date_scaled[j, k, l] + # a spatiotemporally specific intercept
            p_species_date_sq[i]*(date_scaled[j, k, l])^2 # a spatiotemporally specific intercept
          
        }
      }
    }  
  }

  
  # generate initial presence/absence states
  z[,,1] <- rbinom(n=n_sites, size=1, prob=psi1) #
  sum(z[,,1]) / n_sites # true occupancy proportion in year 1
  
  # generate p with heterogeneity
  logit_phi <- array(NA, dim = c(n_species, n_sites, n_years)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        
        logit_phi[i,j,k] = 
          phi_species[i] +
          phi_habitat * habitat_type[j] 
        
      }
    } 
  }

  
  
  # generate presence/absence in subsequent years
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 2:n_years){
        
        # use z as a switch so we are estimating 
        exp_z <- z[i,j,k-1] * ilogit(logit_phi[i,j,k]) + # survival if z=1
          (1 - z[i,j,k-1]) * gamma # or colonization if z=0
        
        # and then assign z stochastically
        # some sites may transition if they are colonized or local extinction occurs
        # but might otherwise retain their state across years
        z[i,j,k] <- rbinom(n = 1, size = 1, prob = exp_z) 
        
      }
    }    
  }

  
  apply(z, 2, sum) / n_sites # true occupancy proportions in each year
  
  apply(z[,1:(n_sites/2),], 3, sum) # true occupancy proportions in each year
  apply(z[,(n_sites/2+1):n_sites,], 3, sum) # true occupancy proportions in each year
  
  # detection / non-detection data
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_visits){
          y[i,j,k,l] <- rbinom(n = 1, size = 1, prob = z[i,j,k]*ilogit(logit_p[i,j,k,l]))
        }
      }
    }
  }

  y ; str(y)
  
  sum((y/n_visits) / sum(z)) # proportion of times detection given presence
  
  # need to fix this if going to use
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
  psi_habitat0 <- numeric(n_years) ; psi_habitat0[1] <- psi1
  psi_habitat1 <- numeric(n_years) ; psi_habitat1[1] <- psi1
  
  # expected for the average species across the two habitat types
    for(k in 2:n_years){ # compute true values of psi
      psi_habitat0[k] <- psi_habitat0[k-1] * ilogit(phi0) + (1 - psi_habitat0[k-1]) * gamma
    }
  
    for(k in 2:n_years){ # compute true values of psi
      psi_habitat1[k] <- psi_habitat1[k-1] * ilogit(phi0 + phi_habitat) + (1 - psi_habitat1[k-1]) * gamma
    }
  
  # actual occupancy rate given the simulated data (for the average species)
  # only works if 1 species # psi_fs_habitat0 <- colSums(z[1:n_species,1:(n_sites/2),]) / (n_sites/2) # psi finite sample
  psi_fs_habitat0 <- apply(z[,1:(n_sites/2),], 3, function(x) mean(x))
  # only works if 1 species # psi_fs_habitat1 <- colSums(z[,(n_sites/2+1):n_sites,]) / (n_sites/2) # psi finite sample
  psi_fs_habitat1 <- apply(z[,(n_sites/2+1):n_sites,], 3, function(x) mean(x))
  
  # compute observed occupancy proportion (for the average species)
  zobs <- apply(y2, c(1,2,3), function(x) max(x, na.rm = TRUE))
  zobs[zobs == "-Inf"] <- NA
  psi_obs_habitat0 <- apply(zobs[,1:(n_sites/2),], 3, sum, na.rm = TRUE) / apply(zobs, 3, function(x) sum(!is.na(x)))
  psi_obs_habitat1 <- apply(zobs[,(n_sites/2+1):n_sites,], 3, sum, na.rm = TRUE) / apply(zobs, 3, function(x) sum(!is.na(x)))
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    psi_habitat0 = psi_habitat0,
    psi_habitat1 = psi_habitat1,
    psi_fs_habitat0 = psi_fs_habitat0,
    psi_fs_habitat1 = psi_fs_habitat1,
    psi_obs_habitat0 = psi_obs_habitat0,
    psi_obs_habitat1 = psi_obs_habitat1,
    psi_eq_habitat0 = psi_eq_habitat0,
    psi_eq_habitat1 = psi_eq_habitat1,
    V = y2, # return detection data after potentially introducing NAs,
    habitat_type = habitat_type,
    date_scaled = date_scaled
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

set.seed(1)
my_simulated_data <- simulate_data(  
  n_species, n_sites, n_years, n_visits,
  psi1, phi0, sigma_phi_species, phi_habitat, gamma, 
  p0, sigma_p_species, p_habitat_type,
  mu_p_species_date, sigma_p_species_date, mu_p_species_date_sq, sigma_p_species_date_sq,
  mean_survey_date, sigma_survey_date,
  create_missing_data, prob_missing
  )

V <- my_simulated_data$V
habitat_type <- my_simulated_data$habitat_type
psi_habitat0 <- my_simulated_data$psi_habitat0
psi_habitat1 <- my_simulated_data$psi_habitat1
psi_fs_habitat0 <- my_simulated_data$psi_fs_habitat0
psi_fs_habitat1 <- my_simulated_data$psi_fs_habitat1
psi_obs_habitat0 <- my_simulated_data$psi_obs_habitat0
psi_obs_habitat1 <- my_simulated_data$psi_obs_habitat1
psi_eq_habitat0 <- my_simulated_data$psi_eq_habitat0
psi_eq_habitat1 <- my_simulated_data$psi_eq_habitat1
date_scaled <- my_simulated_data$date_scaled
species <- seq(1, n_species, by=1)


# Plot trajectories of psi, psi_fs, psi_eq and psi_obs
plot(1:n_years, psi_habitat0, type = 'l', col = "red3", lwd = 3, ylim = c(0,1), xlab = "Year",
     ylab = "Occupancy of average species", frame = F)
lines(psi_habitat1, type = 'l', col = "blue", lwd = 3)
points(1:n_years, psi_fs_habitat0, type = 'b', col = "red3", pch = 16, cex = 2)
points(1:n_years, psi_fs_habitat1, type = 'b', col = "blue", pch = 16, cex = 2)
abline(h = psi_eq_habitat0, lty = 2, col = "red3")
abline(h = psi_eq_habitat1, lty = 2, col = "blue")
lines(1:n_years, psi_obs_habitat0, col = "red", lwd = 3)
lines(1:n_years, psi_obs_habitat1, col = "lightblue", lwd = 3)

legend('topright', c('True psi hab0', 'True psi hab1', 
                     'Finite-sample psi hab0', 'Finite-sample psi hab1', 
                     'Equilibrium psi hab0', 'Equilibrium psi hab1', 
                     'Observed psi hab0', 'Observed psi hab1'),
       col = c('red3', 'blue', 
               'red3', 'blue',
               'red3', 'blue',
               'red', 'lightblue'
               ), pch = c(NA, NA, 16, 16, NA, NA, NA, NA), lty = c(1,1,1,1,2,2,1,1),
       lwd = c(3,3,1,1,2,2,3,3), cex = 1.2, bty = 'n')


## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "species",
               "n_species", "n_sites", "n_years", "n_visits",
               "habitat_type", "date_scaled") 

## Parameters monitored
params <- c("psi1", "phi0", "sigma_phi_species", "phi_habitat", "gamma", 
            "p0", "sigma_p_species", "p_habitat", 
            "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq",
            "T_rep", "T_obs", "P_species")

# MCMC settings
n_iterations <- 800
n_thin <- 1
n_burnin <- 400
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  psi1, phi0, sigma_phi_species, phi_habitat, gamma, 
  p0, sigma_p_species, p_habitat, 
  mu_p_species_date, sigma_p_species_date, mu_p_species_date_sq, sigma_p_species_date_sq,
  NA, NA, NA
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1 = runif(1, 0, 1),
       phi0 = runif(1, -1, 1),
       sigma_phi_species = runif(1, 0, 1),
       phi_habitat = runif(1, -1, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       p_habitat = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

targets <- as.data.frame(cbind(params, parameter_values))


## --------------------------------------------------
### Run model

library(rstan)
stan_model <- "./dynamic_occupancy_model/models/dynocc_model.stan"

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
      pars = c("psi1", "phi0", "sigma_phi_species", "phi_habitat", "gamma", 
               "p0", "sigma_p_species", "p_habitat", 
               "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq"))

saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "phi0", "sigma_phi_species", "phi_habitat", "gamma", "psi1"
))

traceplot(stan_out_sim, pars = c(
  "p0", "sigma_p_species", "p_habitat", "p_date", "p_date_sq"
))


### PPC's

print(stan_out_sim, digits = 3, pars = c("P_species"))

fit_summary <- rstan::summary(stan_out_sim)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[308:457,1]))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)

# Evaluation of fit 
# species 1
plot(list_of_draws[,161], list_of_draws[,11], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 1000),
     xlim = c(0, 1000))

abline(0, 1, lwd = 2, col = "black")

# site 2
plot(list_of_draws[,], list_of_draws[,], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

