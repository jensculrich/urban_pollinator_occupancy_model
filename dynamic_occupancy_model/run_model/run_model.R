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
clade <- my_data$clade
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

mean(species_interaction_metrics$d_supplemented_genus)
sd(species_interaction_metrics$d_supplemented_genus)

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

# final_model.stan is currently the final form
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
stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out_no_specialization_effects.rds")

## --------------------------------------------------
### Model diagnostics (PPCs are in a separate R file)

print(stan_out, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3, 
      pars = c("psi1_0", "sigma_psi1_species",
               "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
               "psi1_interaction_1", "psi1_interaction_2",
               
               "gamma0", "sigma_gamma_species",
               "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
               "gamma_interaction_1", "gamma_interaction_2",
               
               "phi0", "sigma_phi_species",
               "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
               "phi_interaction_1", "phi_interaction_2",
               
               "gamma_year", "phi_year"
      ))

print(stan_out, digits = 3, 
      pars = c( "p0", "sigma_p_species", 
                "p_specialization",
                "p_flower_abundance_any",
                "mu_p_species_date", "sigma_p_species_date", 
                "mu_p_species_date_sq", "sigma_p_species_date_sq",
                "p_year"
                
      ))

print(stan_out, digits = 3, 
      pars = c( "psi1_woody_flowers"
                
      ))

# for continuous model
traceplot(stan_out, pars = c(
  "psi1_0","sigma_psi1_species",
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
  "p0", "sigma_p_species", 
  "p_specialization",
  "p_flower_abundance_any",
  "mu_p_species_date", "sigma_p_species_date", 
  "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))


pairs(stan_out,  
          pars = c(
            "psi1_0","sigma_psi1_species",
            "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
            "psi1_interaction_1", "psi1_interaction_2"
          ))


## --------------------------------------------------
### Gather data for species summary table (appendix s2)

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

phenology_peak_mean_estimate <- fit_summary$summary[577:684,1]
phenology_decay_mean_estimate <- fit_summary$summary[685:792,1]

# retrieve total detections per species
total_detections <- vector(length = n_species)
for(i in 1:n_species){
  total_detections[i] <- sum(V[i,,,])
}
  

appendix_s2 <- as.data.frame(cbind(clade, species_names, total_detections,
            phenology_peak_mean_estimate, phenology_decay_mean_estimate, 
            species_interaction_metrics))
row.names(appendix_s2) <- NULL

write.csv(appendix_s2, "./dynamic_occupancy_model/model_outputs/appendix_s2.csv")

cor(appendix_s2$phenology_peak_mean_estimate, appendix_s2$d_scaled_supplemented_genus)
cor(appendix_s2$phenology_decay_mean_estimate, appendix_s2$d_scaled_supplemented_genus)

library(ggplot2)

summary(m1 <- lm(data = appendix_s2, d_scaled_supplemented_genus ~ phenology_peak_mean_estimate))
p <- ggplot(appendix_s2, aes(phenology_peak_mean_estimate, d_scaled_supplemented_genus)) +
  geom_point() + 
  #geom_smooth(alpha=0.3, method="lm") +
  theme_classic() + 
  ylab("specialization (d')") +
  xlab("species effect on phenological peak \n(greater value == later peak date)") +
  theme(axis.text.x = element_text(size=12),
      axis.title.x = element_text(size=14),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14))
p

summary(m2 <- lm(data = appendix_s2, d_scaled_supplemented_genus ~ phenology_decay_mean_estimate))
q <- ggplot(appendix_s2, aes(phenology_decay_mean_estimate, d_scaled_supplemented_genus)) +
  geom_point() + 
  #geom_smooth(alpha=0.3, method="lm") +
  theme_classic() + 
  ylab("specialization (d')") +
  xlab("species-specific phenological decay \n(greater value == longer flight season)") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14))
q

summary(m3 <- lm(data = appendix_s2, d_scaled_supplemented_genus ~ total_detections))
r <- ggplot(appendix_s2, aes(total_detections, d_scaled_supplemented_genus)) +
  geom_point() + 
  #geom_smooth(alpha=0.3, method="lm") +
  theme_classic() + 
  ylab("specialization (d')") +
  xlab("total detections") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14))
r
