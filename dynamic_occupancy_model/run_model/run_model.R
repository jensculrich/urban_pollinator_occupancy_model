### Run dynamic occupancy model using real data from urban parks pollinator communities
# jcu; started December, 2023

library(rstan)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 1 # >=
my_data <- process_raw_data(min_unique_detections)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_years <- my_data$n_years # number of surveys 
n_visits <- my_data$n_visits
n_years_minus1 <- n_years - 1
species <- seq(1, n_species, by=1)
sites <- seq(1, n_sites, by=1)
years <- seq(1, n_years_minus1, by=1)
date_scaled <- my_data$date_scaled
habitat_type <- my_data$habitat_category

species_names <- my_data$species
site_names <- my_data$sites

## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "species", "sites", "years",
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_visits",
               "habitat_type", "date_scaled") 

## Parameters monitored
params <- c("psi1_0",  "sigma_psi1_species", "mu_psi1_habitat", "sigma_psi1_habitat",
            "gamma0",  "sigma_gamma_species", "mu_gamma_habitat", "sigma_gamma_habitat", 
            "phi0", "sigma_phi_species", "mu_phi_habitat", "sigma_phi_habitat", 
            "p0", "sigma_p_species", "p_habitat", # "sigma_p_site", 
            "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq",
            "species_richness", "avg_species_richness_control", "avg_species_richness_enhanced", 
            "psi_eq_habitat0", "psi_eq_habitat1",
            "T_rep", "T_obs", "P_species")

## Parameters monitored (multinormal model)
params <- c("rho", 
            "sigma_species", "species_intercepts",
            "psi1_0",  "mu_psi1_habitat", "sigma_psi1_habitat", "sigma_psi1_site",
            "gamma0", "mu_gamma_habitat", "sigma_gamma_habitat", "sigma_gamma_species",  "sigma_gamma_site",
            "phi0", "mu_phi_habitat", "sigma_phi_habitat", "sigma_phi_species", "sigma_phi_site",
            "p0", "p_habitat", # "sigma_p_site", 
            "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq",
            "species_richness", "avg_species_richness_control", "avg_species_richness_enhanced", 
            "psi_eq_habitat0", "psi_eq_habitat1",
            "T_rep", "T_obs", "P_species")

# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       #sigma_psi1_species = runif(1, 0, 1),
       mu_psi1_habitat = runif(1, -1, 1),
       sigma_psi1_habitat = runif(1, 0, 1),
       gamma0 = runif(1, -1, 1),
       #sigma_gamma_species = runif(1, 0, 1),
       mu_gamma_habitat = runif(1, -2, -1), # colonization rates are usually low
       sigma_gamma_habitat = runif(1, 0, 1),
       phi0 = runif(1, -1, 1),
       #sigma_phi_species = runif(1, 0, 1),
       mu_phi_habitat = runif(1, 0, 1), # persistence rates are usually greater than 50%
       sigma_phi_habitat = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       #sigma_p_species = runif(1, 0, 1),
       #sigma_p_site = runif(1, 0, 1),
       p_habitat = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

## --------------------------------------------------
### Run model

# stan_model <- "./dynamic_occupancy_model/models/dynocc_model_with_year_effects.stan"
# stan_model <- "./dynamic_occupancy_model/models/dynocc_model.stan"
stan_model <- "./dynamic_occupancy_model/models/dynocc_model_multinormal_w_site_REs.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     control=list(adapt_delta=delta),
                     open_progress = FALSE,
                     cores = n_cores)

saveRDS(stan_out, "./dynamic_occupancy_model/model_outputs/stan_out_multinormal.rds")
stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")


print(stan_out, digits = 3, 
      pars = c("psi1_0",  "sigma_psi1_species", "mu_psi1_habitat", "sigma_psi1_habitat",
               "gamma0",  "sigma_gamma_species", "mu_gamma_habitat", "sigma_gamma_habitat",
               "phi0", "sigma_phi_species", "mu_phi_habitat", "sigma_phi_habitat", 
               "p0", "sigma_p_species", "p_habitat", #  "sigma_p_site",
               "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq"
      ))

print(stan_out, digits = 3, 
      pars = c("species_richness", "avg_species_richness_control", "avg_species_richness_enhanced"
      ))

print(stan_out, digits = 3, 
      pars = c("psi_eq_habitat0", "psi_eq_habitat1"
      ))

print(stan_out, digits = 3, 
      pars = c("gamma_year", "phi_year"
      ))

traceplot(stan_out, pars = c(
  "psi1_0",  #"sigma_psi1_species", 
  "mu_psi1_habitat", "sigma_psi1_habitat", "sigma_psi1_site",
  "gamma0",  "sigma_gamma_species", 
  "mu_gamma_habitat", "sigma_gamma_habitat", "sigma_gamma_site",
  "phi0", "sigma_phi_species", 
  "mu_phi_habitat", "sigma_phi_habitat", "sigma_phi_site"
))


traceplot(stan_out, pars = c(
  "rho", 
  "sigma_species"
))

traceplot(stan_out, pars = c(
  "gamma_year",  "phi_year"
))

traceplot(stan_out, pars = c(
  "p0", #"sigma_p_species", 
  #"sigma_p_site", 
  "p_habitat", 
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))

traceplot(stan_out,  
      pars = c("avg_species_richness_control", "avg_species_richness_enhanced"
      ))

print(stan_out, digits = 3, pars = c("P_species"))

# get an "average" P value
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[478:575,1]))
