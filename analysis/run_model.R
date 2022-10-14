### Run occupancy model using real data from urban parks pollinator communities
# jcu; started oct 13, 2022

library(rstan)

source("./analysis/prep_data.R")
my_data <- process_raw_data()

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 
n_visits <- my_data$n_visits

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

date_scaled <- my_data$date_scaled
habitat_category <- my_data$habitat_category

stan_data <- c("V", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites", 
               "habitat_category", "date_scaled") 

# Parameters monitored
params <- c("mu_psi_0",
            "sigma_psi_species",
            "mu_psi_species_habitat",
            "sigma_psi_species_habitat",
            "mu_p_0",
            "sigma_p_species",
            "mu_p_species_habitat",
            "sigma_p_species_habitat",
            "mu_p_species_date",
            "sigma_p_species_date",
            "mu_p_species_date_sq",
            "sigma_p_species_date_sq"
)


# MCMC settings
n_iterations <- 200
n_thin <- 1
n_burnin <- 100
n_chains <- 3
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       mu_psi_species_habitat = runif(1, -1, 1),
       sigma_psi_species_habitat = runif(1, 0, 1),
       mu_p_0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       mu_p_species_habitat = runif(1, -1, 1),
       sigma_p_species_habitat = runif(1, 0, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

#targets <- as.data.frame(cbind(param_names, parameters))

## --------------------------------------------------
### Run model

stan_model <- "./models/model0.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

print(stan_out, digits = 3)
saveRDS(stan_out, "./analysis/stan_out.rds")
stan_out <- readRDS("./analysis/stan_out.rds")