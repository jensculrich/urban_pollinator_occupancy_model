
## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_species <- 120 # number of species
n_sites <- 60 # number of sites (must be an even number for simulation code)
n_years <- 7 # number of years
n_years_minus1 <- n_years - 1
n_visits <- 6 # number of surveys per year

# set parameter values
psi1_0 <- 0 # prob of initial occupancy
sigma_psi1_species <- 2 # prob of initial occupancy
mu_psi1_habitat <- 0.5 # habitat effect
sigma_psi1_habitat <- 1 # species-specific variation

gamma0 <- -2 # colonization probability
sigma_gamma_species <- 0.5 # species-specific variation
mu_gamma_habitat <- 0.25 # habitat effect
sigma_gamma_habitat <- 0.75 # species-specific variation
# add year effects (logit-scale) # make min and max == 0 if not including this in model
# min must be lower than max
gamma_min = -0.5
gamma_max = 0.5

phi0 <- 1.25 # persistence probability
sigma_phi_species <- 1 # species-specific variation
mu_phi_habitat <- 1.5 # habitat effect
sigma_phi_habitat <- 1 # species-specific variation
# add year effects (logit-scale) # make min and max == 0 if not including this in model
# min must be lower than max
phi_min = -1
phi_max = 1

p0 <- 0 # probability of detection (logit scaled)
sigma_p_species <- 2 # species-specific variation
sigma_p_site <- 0.5 # site-specific variation
# site specific random effect
p_habitat <- -0.5 # increase in detection rate moving from one habitat type to the other (logit scaled)
mu_p_species_date <- 0
sigma_p_species_date <- 1
mu_p_species_date_sq <- -0.5  
sigma_p_species_date_sq <- 1

mean_survey_date <- 180
sigma_survey_date <- 40

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- FALSE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

# test generating some correlated data
rho <- .8
# Correlation matrix
(rho_u <- matrix(c(1, rho, rho, 1), ncol = 2))

# Cholesky factor:
# (Transpose it so that it looks the same as in Stan)
(L_u <- t(chol(rho_u)))

# Verify that we recover rho_u,
# Recall that %*% indicates matrix multiplication
L_u %*% t(L_u)

# Generate uncorrelated z from a standard normal distribution assuming 50 species
n_species <- 50
delta0_u1 = 0.5
delta1_u1 = 0
delta0_u2 = 0
delta1_u2 = 1

d <- runif(n_species, -2, 2)

mu_z_u1 = vector(length = n_species)

for(i in 1:n_species){
  mu_z_u1[i] = delta0_u1 + delta1_u1*d[i]
}

z_u1 = vector(length = n_species)

for(i in 1:n_species){
  z_u1[i] <- rnorm(n=1, mean=mu_z_u1[i], sd=1)
}

mu_z_u2 = vector(length = n_species)

for(i in 1:n_species){
  mu_z_u2[i] = delta0_u2 + delta1_u2*d[i]
}

z_u2 = vector(length = n_species)

for(i in 1:n_species){
  z_u2[i] <- rnorm(n=1, mean=mu_z_u2[i], sd=1)
}

#(z_u1 <- rnorm(n_species, 0, 1))
#(z_u2 <- rnorm(n_species, 0, 1))

# matrix z_u
(z_u <- matrix(c(z_u1, z_u2), ncol = n_species, byrow = TRUE))

# Then, generate correlated parameters by pre-multiplying the z_u matrix with the L_u
L_u %*% z_u

# Use the following diagonal matrix to scale the z_u.
tau_u1 <- 0.2
tau_u2 <- 0.1
(diag_matrix_tau <- diag(c(tau_u1,  tau_u2)))

# Finally, generate the adjustments for each subject u:
(u <- diag_matrix_tau %*% L_u %*% z_u)

# The rows are correlated ~.8
cor(u[1, ], u[2, ])

# The variance components can be recovered as well:
sd(u[1, ])
sd(u[2, ])

plot(u[1, ], d)
cor(u[1, ], d)

plot(u[2, ], d)
cor(u[2, ], d)

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
    n_species, n_sites, n_years, n_visits,
    psi1_0, sigma_psi1_species, mu_psi1_habitat, sigma_psi1_habitat,
    gamma0, sigma_gamma_species, mu_gamma_habitat, sigma_gamma_habitat, 
    phi0, sigma_phi_species, mu_phi_habitat, sigma_phi_habitat, 
    p0, sigma_p_species, sigma_p_site, p_habitat_type,
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
  psi1_0 <- psi1_0 # prob of initial occupancy
  sigma_psi1_species <- sigma_psi1_species # species-specific variation
  mu_psi1_habitat <- mu_psi1_habitat # effect of habitat type on occurrence
  sigma_psi1_habitat <- sigma_psi1_habitat
  
  gamma0 <- gamma0 # colonization probability
  sigma_gamma_species <- sigma_gamma_species # species-specific variation
  mu_gamma_habitat <- mu_gamma_habitat # effect of habitat type on colonization
  sigma_gamma_habitat <- sigma_gamma_habitat # effect of habitat type on colonization
  gamma_year <- vector(length = n_years_minus1)
  for(k in 1:length(gamma_year)){
    gamma_year[k] <- runif(1, min=gamma_min, max=gamma_max)
  }
  
  phi0 <- phi0 # persistence probability
  sigma_phi_species <- sigma_phi_species # species-specific variation
  mu_phi_habitat <- mu_phi_habitat # effect of habitat type on persistence
  sigma_phi_habitat <- sigma_phi_habitat # effect of habitat type on persistence
  phi_year <- vector(length = n_years_minus1)
  for(k in 1:length(phi_year)){
    phi_year[k] <- runif(1, min=phi_min, max=phi_max)
  }
  
  # p <- p # probability of detection
  # equilibrium occupancy rate (for the average species)
  (psi_eq_habitat0 <- ilogit(gamma0) / (ilogit(gamma0)+(1-ilogit(phi0)))) # equilibrium occupancy rate 
  (psi_eq_habitat1 <- ilogit(gamma0 + mu_gamma_habitat) / (ilogit(gamma0 + mu_gamma_habitat)+(1-ilogit(phi0 + mu_phi_habitat)))) # equilibrium occupancy rate
  
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
  psi1_species <- rnorm(n=n_species, mean=psi1_0, sd=sigma_psi1_species)
  psi1_habitat <- rnorm(n=n_species, mean=mu_psi1_habitat, sd=sigma_psi1_habitat)
  
  gamma_species <- rnorm(n=n_species, mean=gamma0, sd=sigma_gamma_species)
  gamma_habitat <- rnorm(n=n_species, mean=mu_gamma_habitat, sd=sigma_gamma_habitat)
  
  phi_species <- rnorm(n=n_species, mean=phi0, sd=sigma_phi_species)
  phi_habitat <- rnorm(n=n_species, mean=mu_phi_habitat, sd=sigma_phi_habitat)
  
  p_species <- rnorm(n=n_species, mean=p0, sd=sigma_p_species)

  p_site <- rnorm(n=n_sites, mean=0, sd=sigma_p_site)
  
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
            p_site[j] +
            p_habitat * habitat_type[j] +
            p_species_date[i]*date_scaled[j, k, l] + # a spatiotemporally specific intercept
            p_species_date_sq[i]*(date_scaled[j, k, l])^2 # a spatiotemporally specific intercept
          
        }
      }
    }  
  }

  
  # generate ecological expected values with heterogeneity
  logit_psi1 <- array(NA, dim = c(n_species, n_sites)) 
  logit_gamma <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  logit_phi <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years_minus1){
        
        logit_psi1[i,j] = 
          psi1_species[i] +
          psi1_habitat[i] * habitat_type[j] 
        
        logit_gamma[i,j,k] = 
          gamma_species[i] +
          gamma_habitat[i] * habitat_type[j] +
          gamma_year[k]
        
        logit_phi[i,j,k] = 
          phi_species[i] +
          phi_habitat[i] * habitat_type[j] +
          phi_year[k]
        
      }
    } 
  }

  # generate initial presence/absence states
  for(i in 1:n_species){
    for(j in 1:n_sites){
      z[i,j,1] <- rbinom(n=1, size=1, prob=ilogit(logit_psi1[i,j])) 
    }
  }

  sum(z[,,1]) / n_sites # true occupancy proportion in year 1
  
  # generate presence/absence in subsequent years
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 2:n_years){
        
        # use z as a switch so we are estimating 
        exp_z <- z[i,j,k-1] * ilogit(logit_phi[i,j,k-1]) + # survival if z=1
          (1 - z[i,j,k-1]) * ilogit(logit_gamma[i,j,k-1]) # or colonization if z=0
        
        # and then assign z stochastically
        # some sites may transition if they are colonized or local extinction occurs
        # but might otherwise retain their state across years
        z[i,j,k] <- rbinom(n = 1, size = 1, prob = exp_z) 
        
      }
    }    
  }

  apply(z, 3, sum) / n_sites # true occupancy proportions in each year
  
  apply(z[,1:(n_sites/2),], 3, sum) / (n_sites / 2) # true occupancy proportions in each year
  apply(z[,(n_sites/2+1):n_sites,], 3, sum) / (n_sites / 2) # true occupancy proportions in each year
  
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
  
  # which species never occurred
  species_not_occurring <- length(which(apply(z, 1, sum) == 0))
  
  # which species were never detected (even though they occurred)
  species_not_observed <- length(which(apply(y2, 1, sum) == 0))
  
  table(nvisits <- apply(y2, c(1,3), function(x) sum(!is.na(x))))
  # _ sites with no visits, _ with 1 visit, and _ with 2 visits through 7 visits
  
  # compute true expected and realized occupancy (psi and psi_fs)
  psi_habitat0 <- numeric(n_years) ; psi_habitat0[1] <- ilogit(psi1_0)
  psi_habitat1 <- numeric(n_years) ; psi_habitat1[1] <- ilogit(psi1_0 + mu_psi1_habitat)
  
  # expected for the average species across the two habitat types
    for(k in 2:n_years){ # compute true values of psi
      psi_habitat0[k] <- psi_habitat0[k-1] * ilogit(phi0) + (1 - psi_habitat0[k-1]) * ilogit(gamma0)
    }
  
    for(k in 2:n_years){ # compute true values of psi
      psi_habitat1[k] <- psi_habitat1[k-1] * ilogit(phi0 + mu_phi_habitat) + (1 - psi_habitat1[k-1]) * ilogit(gamma0 + mu_gamma_habitat)
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
  
  print(paste0("*** ", species_not_occurring, " / ", n_species, " species never occurred ***"))
  print(paste0("*** ", species_not_observed, " / ", n_species, " species never detected ***"))
  
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
    date_scaled = date_scaled,
    gamma_year = gamma_year,
    phi_year = phi_year
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

set.seed(1)
my_simulated_data <- simulate_data(  
  n_species, n_sites, n_years, n_visits,
  psi1_0, sigma_psi1_species, mu_psi1_habitat, sigma_psi1_habitat,
  gamma0, sigma_gamma_species, mu_gamma_habitat, sigma_gamma_habitat, 
  phi0, sigma_phi_species, mu_phi_habitat, sigma_phi_habitat, 
  p0, sigma_p_species, sigma_p_site, p_habitat_type,
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
sites <- seq(1, n_sites, by=1)
years <- seq(1, n_years_minus1, by=1)

gamma_year <- my_simulated_data$gamma_year
phi_year <- my_simulated_data$phi_year

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

stan_data <- c("V", "species", "sites", "years",
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_visits",
               "habitat_type", "date_scaled") 

## Parameters monitored
params <- c("psi1_0",  "sigma_psi1_species", "mu_psi1_habitat", "sigma_psi1_habitat",
            "gamma0",  "sigma_gamma_species", "mu_gamma_habitat", "sigma_gamma_habitat", "gamma_year",
            "phi0", "sigma_phi_species", "mu_phi_habitat", "sigma_phi_habitat", "phi_year",
            "p0", "sigma_p_species", "sigma_p_site", "p_habitat", 
            "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq",
            "T_rep", "T_obs", "P_species")

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  psi1_0, sigma_psi1_species, mu_psi1_habitat, sigma_psi1_habitat,
  gamma0, sigma_gamma_species, mu_gamma_habitat, sigma_gamma_habitat, NA,
  phi0, sigma_phi_species, mu_phi_habitat, sigma_phi_habitat, NA,
  p0, sigma_p_species, sigma_p_site, p_habitat, 
  mu_p_species_date, sigma_p_species_date, mu_p_species_date_sq, sigma_p_species_date_sq,
  NA, NA, NA
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0, 1),
       mu_psi1_habitat = runif(1, -1, 1),
       sigma_psi1_habitat = runif(1, 0, 1),
       gamma0 = runif(1, -1, 1),
       sigma_gamma_species = runif(1, 0, 1),
       mu_gamma_habitat = runif(1, -2, -1), # colonization rates are usually low
       sigma_gamma_habitat = runif(1, 0, 1),
       phi0 = runif(1, -1, 1),
       sigma_phi_species = runif(1, 0, 1),
       mu_phi_habitat = runif(1, 0, 1), # persistence rates are usually greater than 50%
       sigma_phi_habitat = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
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
stan_model <- "./dynamic_occupancy_model/models/dynocc_model_with_year_effects.stan"

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
      pars = c("psi1_0",  "sigma_psi1_species", "mu_psi1_habitat", "sigma_psi1_habitat",
               "gamma0",  "sigma_gamma_species", "mu_gamma_habitat", "sigma_gamma_habitat",
               "phi0", "sigma_phi_species", "mu_phi_habitat", "sigma_phi_habitat", 
               "p0", "sigma_p_species", "sigma_p_site", "p_habitat", 
               "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq"
               ))

print(stan_out_sim, digits = 3, 
      pars = c("gamma_year", "phi_year"
      ))

print(gamma_year)
print(phi_year)

saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "psi1_0",  "sigma_psi1_species", "mu_psi1_habitat", "sigma_psi1_habitat",
  "gamma0",  "sigma_gamma_species", "mu_gamma_habitat", "sigma_gamma_habitat",
  "phi0", "sigma_phi_species", "mu_phi_habitat", "sigma_phi_habitat"
))

traceplot(stan_out_sim, pars = c(
  "gamma_year",  "phi_year"
))

traceplot(stan_out_sim, pars = c(
  "p0", "sigma_p_species", "sigma_p_site", "p_habitat", 
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))


### PPC's

print(stan_out_sim, digits = 3, pars = c("P_species"))

fit_summary <- rstan::summary(stan_out_sim)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[:,1]))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)

# Evaluation of fit 
# species 1
plot(list_of_draws[,127], list_of_draws[,27], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))

abline(0, 1, lwd = 2, col = "black")

# discrepancy in the actual data too low (or discrepancy in replicate data too high)?

# site 2
plot(list_of_draws[,], list_of_draws[,], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

