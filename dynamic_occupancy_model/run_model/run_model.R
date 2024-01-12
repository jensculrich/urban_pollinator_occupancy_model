### Run dynamic occupancy model using real data from urban parks pollinator communities
# jcu; started December, 2023

library(rstan)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 3 # >=
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
site_year_visit_count <- my_data$site_year_visit_count
date_scaled <- my_data$date_scaled
habitat_type <- my_data$habitat_category
species_interaction_metrics <- my_data$species_interaction_metrics
d <- species_interaction_metrics$d_scaled
degree <- species_interaction_metrics$degree_scaled
species_names <- my_data$species
site_names <- my_data$sites
herbaceous_flowers_scaled <- my_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_data$woody_flowers_scaled
herbaceous_flowers_by_survey <- my_data$herbaceous_flowers_by_survey
woody_flowers_by_survey <- my_data$woody_flowers_by_survey
flowers_any_by_survey <- my_data$flowers_any_by_survey

## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "species", "sites", "years",
               "n_species", "n_sites", "n_years", "n_years_minus1", "site_year_visit_count", "n_visits",
               "habitat_type", "date_scaled", "d", "degree", "herbaceous_flowers_scaled", "woody_flowers_scaled", "flowers_any_by_survey") 

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
params <- c(
  "L_psi1_species", "sigma_psi1_species",
  "delta0_psi1_herbaceous", "delta1_psi1_herbaceous", 
  "delta0_psi1_woody", "delta1_psi1_woody", 
  "psi1_0",
  
  "L_gamma_species", "sigma_gamma_species",
  "delta0_gamma_herbaceous", "delta1_gamma_herbaceous", 
  "delta0_gamma_woody", "delta1_gamma_woody", 
  "gamma0",
  
  "phi0", "delta1_phi0", "sigma_phi_species",
  "delta0_phi_herbaceous", "delta1_phi_herbaceous", "sigma_phi_herbaceous",
  "delta0_phi_woody", "delta1_phi_woody", "sigma_phi_woody",
  
  "p0", "delta1_p0", "sigma_p_species", 
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq", "p_flower_abundance_any", 
  "species_richness", "avg_species_richness_control", "avg_species_richness_enhanced", 
  #"turnover_control", "turnover_enhanced",
  #"psi_eq_habitat0", "psi_eq_habitat1",
  "T_rep", "T_obs", "P_species")

## Parameters monitored (multinormal model with species effects and continuous flower plant data)
params <- c(
            "psi1_0", "sigma_psi1_species",
            "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
            "psi1_interaction_1", "psi1_interaction_2",
  
            "gamma0", "sigma_gamma_species",
            "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
            "gamma_interaction_1", "gamma_interaction_2",
            
            "phi0", "sigma_phi_species",
            "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
            "phi_interaction_1", "phi_interaction_2",
            
            "p0", "delta1_p0", "sigma_p_species", 
            "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq", "p_flower_abundance_any", 
            "species_richness", "avg_species_richness_control", "avg_species_richness_enhanced", 
            #"turnover_control", "turnover_enhanced",
            #"psi_eq_habitat0", "psi_eq_habitat1",
            "T_rep", "T_obs", "P_species")

# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 4
n_cores <- n_chains
delta = 0.9

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample

inits <- lapply(1:n_chains, function(i)
  
  list(
       psi1_0 = runif(1, 0, 1),
       delta0_psi1_herbaceous = runif(1, -1, 1),
       delta1_psi1_herbaceous = runif(1, -1, 1),
       delta0_psi1_woody = runif(1, -1, 1), 
       delta1_psi1_woody = runif(1, -1, 1),
       gamma0 = runif(1, -1, 0),
       delta0_gamma_herbaceous = runif(1, -1, 1),
       delta1_gamma_herbaceous = runif(1, -1, 1),
       delta0_gamma_woody = runif(1, -1, 1), 
       delta1_gamma_woody = runif(1, -1, 1),
       phi0 = runif(1, 1, 3),
       delta0_phi_herbaceous = runif(1, -1, 1),
       delta1_phi_herbaceous = runif(1, -1, 1),
       delta0_phi_woody = runif(1, -1, 1), 
       delta1_phi_woody = runif(1, -1, 1),
       p0 = runif(1, -1, 0),
       sigma_p_species = runif(1, 0, 1),
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
stan_model <- "./dynamic_occupancy_model/models/dynocc_model_4.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                     data = stan_data, 
                     #init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     control=list(adapt_delta=delta),
                     open_progress = FALSE,
                     cores = n_cores)

saveRDS(stan_out, "./dynamic_occupancy_model/model_outputs/stan_out2.rds")
stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")


print(stan_out, digits = 3, 
      pars = c("psi1_0", "delta1_psi1_0", "sigma_psi1_species",   
               "delta0_psi1_herbaceous", "delta1_psi1_herbaceous",
               "delta0_psi1_woody", "delta1_psi1_woody",
               "sigma_psi1_herbaceous", "sigma_psi1_woody", 
               "gamma0", "delta1_gamma0", "sigma_gamma_species",
               "delta0_gamma_herbaceous", "delta1_gamma_herbaceous",
               "delta0_gamma_woody", "delta1_gamma_woody",
               "sigma_gamma_herbaceous", "sigma_gamma_woody", 
               "phi0", "delta1_phi0", "sigma_phi_species",
               "delta0_phi_herbaceous", "delta1_phi_herbaceous",
               "delta0_phi_woody", "delta1_phi_woody",
               "sigma_phi_herbaceous", "sigma_phi_woody"
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

print(stan_out, digits = 3, 
      pars = c("delta0_phi0", "delta1_phi0",
               "delta0_phi_habitat", "delta1_phi_habitat"
      ))

traceplot(stan_out, pars = c(
  "psi1_0",  
  "delta0_psi1_habitat", "delta1_psi1_habitat",
  "sigma_psi1_habitat", 
  "gamma0", "delta1_gamma0",
  "delta0_gamma_habitat", "delta1_gamma_habitat",
  "epsilon0_gamma_habitat", "epsilon1_gamma_habitat",
  "sigma_gamma_species", "sigma_gamma_habitat", 
  "phi0", "delta1_phi0",
  "delta0_phi_habitat", "delta1_phi_habitat",
  "epsilon0_phi_habitat", "epsilon1_phi_habitat",
  "sigma_phi_habitat" , 
  "sigma_phi_species"
))

# for continuous model
traceplot(stan_out, pars = c(
  "psi1_0", "delta1_psi1_0", "sigma_psi1_species",   
  "delta0_psi1_herbaceous", "delta1_psi1_herbaceous",
  "delta0_psi1_woody", "delta1_psi1_woody",
  "sigma_psi1_herbaceous", "sigma_psi1_woody", 
  "gamma0", "delta1_gamma0", "sigma_gamma_species",
  "delta0_gamma_herbaceous", "delta1_gamma_herbaceous",
  "delta0_gamma_woody", "delta1_gamma_woody",
  "sigma_gamma_herbaceous", "sigma_gamma_woody", 
  "phi0", "delta1_phi0", "sigma_phi_species",
  "delta0_phi_herbaceous", "delta1_phi_herbaceous",
  "delta0_phi_woody", "delta1_phi_woody",
  "sigma_phi_herbaceous", "sigma_phi_woody"
))

traceplot(stan_out, pars = c(
  "sigma_psi1_species", "L_psi1_species",
  "gamma0", "delta1_gamma0", "sigma_gamma_species",
  "delta0_gamma_herbaceous", "delta1_gamma_herbaceous",
  "delta0_gamma_woody", "delta1_gamma_woody",
  "sigma_gamma_herbaceous", "sigma_gamma_woody",
  "phi0",  "delta1_phi0", "sigma_phi_species",
  "delta0_phi_herbaceous", "delta1_phi_herbaceous",
  "delta0_phi_woody", "delta1_phi_woody",
  "sigma_phi_herbaceous", "sigma_phi_woody"
  
))

# for continuous model
traceplot(stan_out, pars = c(
  "delta0_psi1_herbaceous", "delta1_psi1_herbaceous",
  "delta0_psi1_woody", "delta1_psi1_woody", "psi1_0"
))

# for continuous model
traceplot(stan_out, pars = c(
  "sigma_gamma_species", "L_gamma_species",
  "gamma0",
  "delta0_gamma_herbaceous", "delta1_gamma_herbaceous",
  "delta0_gamma_woody", "delta1_gamma_woody"
))

# for continuous model
traceplot(stan_out, pars = c(
  "psi1_0", "sigma_psi1_species",
  "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
  "psi1_interaction_1", "psi1_interaction_2",
  
  "gamma0", "sigma_gamma_species",
  "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
  "gamma_interaction_1", "gamma_interaction_2",
  
  "phi0", "sigma_phi_species",
  "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
  "phi_interaction_1", "phi_interaction_2"
))

traceplot(stan_out, pars = c(
  "rho", 
  "sigma_species"
))

traceplot(stan_out, pars = c(
  "gamma_year",  "phi_year"
))

traceplot(stan_out, pars = c(
  "p0", "delta1_p0", "sigma_p_species", 
  #"sigma_p_site", 
  #"p_habitat", 
  "p_flower_abundance_any",
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))

traceplot(stan_out,  
      pars = c("avg_species_richness_control", "avg_species_richness_enhanced"
      ))

traceplot(stan_out,  
          pars = c("species_richness"
          ))

print(stan_out, digits = 3, pars = c("P_species"))

# get an "average" P value
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[236:308,1]))
