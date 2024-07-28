### Run dynamic occupancy model using real data from urban parks pollinator communities
# jcu; started December, 2023

library(rstan)

min_unique_detections = 1 
filter_nonnative_woody = FALSE
supplement_interactions = TRUE
cluster_interactions_at_genus = TRUE
source("./dynamic_occupancy_model/run_model/prep_data.R")
my_data <- process_raw_data(min_unique_detections, filter_nonnative_woody)

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
years_full <- seq(1, n_years, by=1)
site_year_visit_count <- my_data$site_year_visit_count
date_scaled <- my_data$date_scaled
date_adjusted_scaled <- my_data$date_adjusted_scaled
habitat_type <- my_data$habitat_category
species_interaction_metrics <- my_data$species_interaction_metrics
# use interaction data from our study only (FALSE)? or supplement with external observations (TRUE)?
if(supplement_interactions == TRUE){
  d <- species_interaction_metrics$d_scaled_supplemented
  degree <- species_interaction_metrics$degree_scaled_supplemented
} else {
  d <- species_interaction_metrics$d_scaled
  degree <- species_interaction_metrics$degree_scaled
}
# replace with genus level specialization?
if(cluster_interactions_at_genus == TRUE){
  d <- species_interaction_metrics$d_scaled_supplemented_genus
  degree <- species_interaction_metrics$degree_scaled_supplemented_genus
}
species_names <- my_data$species
site_names <- my_data$sites
herbaceous_flowers_scaled <- my_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_data$woody_flowers_scaled
herbaceous_diversity_scaled <- my_data$herbaceous_diversity_scaled
woody_diversity_scaled <- my_data$woody_diversity_scaled
woody_diversity_scaled_avg <- cbind(rowMeans(woody_diversity_scaled), rowMeans(woody_diversity_scaled), rowMeans(woody_diversity_scaled))
herbaceous_flowers_by_survey <- my_data$herbaceous_flowers_by_survey
woody_flowers_by_survey <- my_data$woody_flowers_by_survey
flowers_any_by_survey <- my_data$flowers_any_by_survey
woody_flowers_scaled_all_years <- my_data$woody_flowers_scaled_all_years
woody_flowers_scaled_all_years_repped <- cbind(woody_flowers_scaled_all_years, woody_flowers_scaled_all_years, woody_flowers_scaled_all_years)
preestablished <- my_data$preestablished

# write.csv(as.data.frame(species_interaction_metrics), "./data/species_information.csv")

## --------------------------------------------------
### Prep data and tweak model options

# use adjusted dates?
#date_scaled <- date_adjusted_scaled
# use the diversity instead of abundance data
#herbaceous_flowers_scaled <- herbaceous_diversity_scaled
#woody_flowers_scaled <- woody_diversity_scaled_avg

stan_data <- c("V", "species", "sites", "years", "years_full", 
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_visits",
               "habitat_type", "date_scaled", "d", "degree",
               "herbaceous_flowers_scaled", 
               "woody_flowers_scaled", 
               "flowers_any_by_survey"
) 

## Parameters monitored 
params <- c(#"L_species", "sigma_species",
  
  "psi1_0", 
  "sigma_psi1_species",
  "psi1_herbaceous_flowers", "psi1_woody_flowers", 
  "psi1_specialization",
  "psi1_interaction_1", "psi1_interaction_2",

  "gamma0", 
  "sigma_gamma_species",
  "gamma_herbaceous_flowers", "gamma_woody_flowers", 
  "gamma_specialization",
  "gamma_interaction_1", "gamma_interaction_2", 

  "phi0", 
  "sigma_phi_species",
  "phi_herbaceous_flowers", "phi_woody_flowers", 
  "phi_specialization",
  "phi_interaction_1", "phi_interaction_2",

  "p0", 
  "sigma_p_species", 
  "p_specialization",
  "mu_p_species_date", "sigma_p_species_date", 
  "mu_p_species_date_sq", "sigma_p_species_date_sq", 
  "p_flower_abundance_any", 
  
  "gamma_year",
  "phi_year",
  "p_year",
  
  "W_species_rep",
  "psi1_species", "gamma_species", "phi_species", 
  "p_species", "p_date", "p_date_sq")

# MCMC settings
n_iterations <- 4000
n_thin <- 1
n_burnin <- n_iterations/2
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 1, 2),
       psi1_herbaceous_flowers = runif(1, -1, 1),
       psi1_woody_flowers = runif(1, -1, 1),
       psi1_specialization = runif(1, -1, 1),
       psi1_interaction_1 = runif(1, -1, 1),
       psi1_interaction_2 = runif(1, -1, 1),
       
       gamma0 = runif(1, -1, 0),
       sigma_gamma_species = runif(1, 1, 2),
       gamma_herbaceous_flowers = runif(1, -1, 1),
       gamma_woody_flowers = runif(1, -1, 1),
       gamma_specialization = runif(1, -1, 1),
       gamma_interaction_1 = runif(1, -1, 1),
       gamma_interaction_2 = runif(1, -1, 1),
       
       phi0 = runif(1, 1, 2),
       sigma_phi_species = runif(1, 1, 2),
       phi_herbaceous_flowers = runif(1, -1, 1),
       phi_woody_flowers = runif(1, -1, 1),
       phi_specialization = runif(1, -1, 1),
       phi_interaction_1 = runif(1, -1, 1),
       phi_interaction_2 = runif(1, -1, 1),
       
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 1, 2),
       p_degree = runif(1, -1, 1),
       p_flower_abundance_any = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

## --------------------------------------------------
### Run model
# dynocc7 is currently the final form
stan_model <- "./dynamic_occupancy_model/models/final_model.stan"

## Call Stan from R
set.seed(1)
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

saveRDS(stan_out, "./dynamic_occupancy_model/model_outputs/stan_out.rds")
stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")


print(stan_out, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), digits = 3, 
      pars = c("psi1_0", #sigma_psi1_species",
               "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
               "psi1_interaction_1", "psi1_interaction_2",
               
               "gamma0", "sigma_gamma_species",
               "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
               "gamma_interaction_1", "gamma_interaction_2",
               
               "phi0", "sigma_phi_species",
               "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
               "phi_interaction_1", "phi_interaction_2"
      ))

print(stan_out, digits = 3, 
      pars = c( "p0", "sigma_p_species", 
                "p_specialization",
                "p_flower_abundance_any",
                "mu_p_species_date", "sigma_p_species_date", 
                "mu_p_species_date_sq", "sigma_p_species_date_sq" 
                
      ))

print(stan_out, digits = 3, 
      pars = c("avg_species_richness_control", "avg_species_richness_enhanced", 
               "increase_richness_enhanced"
      ))

print(stan_out,  
          pars = c("p_year"
          ))

# for continuous model
traceplot(stan_out, pars = c(
  "psi1_0","sigma_psi1_species",
  "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
  "psi1_interaction_1", "psi1_interaction_2",
  
  #"psi1_herbaceous_diversity", "psi1_interaction_3",
  
  "gamma0", "sigma_gamma_species",
  "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
  "gamma_interaction_1", "gamma_interaction_2",
  
  #"gamma_herbaceous_diversity", "gamma_interaction_3",
  
  "phi0", "sigma_phi_species",
  "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
  "phi_interaction_1", "phi_interaction_2"
  
  #"phi_herbaceous_diversity", "phi_interaction_3"
))

traceplot(stan_out, pars = c(
  "p0", "sigma_p_species", 
  "p_specialization",
  "p_flower_abundance_any",
  "mu_p_species_date", "sigma_p_species_date", 
  "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))

traceplot(stan_out,  
          pars = c("L_species", "sigma_species"
          ))

traceplot(stan_out,  
      pars = c("avg_species_richness_control", "avg_species_richness_enhanced", 
               "increase_richness_enhanced"
))

traceplot(stan_out,  
          pars = c("gamma_year", 
                   "phi_year",
                   "p_year"
                   ))

pairs(stan_out,  
          pars = c(
            "psi1_0","sigma_psi1_species",
            "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
            "psi1_interaction_1", "psi1_interaction_2"
          ))

print(stan_out, digits = 3, pars = c("P_species"))

