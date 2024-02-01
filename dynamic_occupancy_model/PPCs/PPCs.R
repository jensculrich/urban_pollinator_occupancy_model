# how many detections would you expect per species (W)?
# we will compare this to how many detections per species simulated in the 
# generated quantities block of our model (W_rep) for a visual PPC

library(tidyverse)

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

# Summarize V[i,j,k,l] by species -> to get W[i]
W_species = vector(length = n_species)

for(i in 1:n_species){
  W_species[i] = sum(V[i,,,])
}


# for simulated data
W_df <- as.data.frame(cbind(species, W_species)) %>%
  mutate(W_species = as.numeric(W_species))

# for real data
W_df <- as.data.frame(cbind(species_names, W_species)) %>%
  mutate(W_species = as.numeric(W_species))

# get W distributions from model
#stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")
#fit_summary <- rstan::summary(stan_out_sim)
fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest


## --------------------------------------------------
# 

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

start = 1 # which species to start at (hard to see them all at once)
# start at 1, 37, and 73 is pretty good for visualization
n = 36 # how many species to plot (36 is a good number to look at the species in 3 slices)

stan_fit_first_W <- 91 # this changes depending on how many params you tracked

df_estimates <- data.frame(X = numeric(), 
                           Y = numeric(), 
                           lower_95 = numeric(),
                           upper_95 = numeric(),
                           lower_50 = numeric(),
                           upper_50 = numeric()
                           ) 

for(i in 1:n){
  
  row <- c((i + start - 1), 
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),1],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),4],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),8],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),5],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),7])
  
  df_estimates[i,] <- row
  
}

labels=as.vector(c(W_df[start:(start + n - 1),1]))
ylims = c(0,(max(df_estimates$upper_95)+5))
end_point  = 0.5 + nrow(df_estimates) + nrow(df_estimates) - 1 #

par(mar = c(9,4,1,2))
plot(1, type="n",
     xlim=c(start, n + start - 1 + 0.5), 
     xlab="",
     xaxt = "n",
     ylim=ylims, 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Real Detections vs. Model Expectations of Detections")

#axis(1, at=start:(start+n), labels=labels, las = 2, cex.axis=.75)
text(seq(start, start + n - 1, by = 1), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = labels, cex = 1)

for(i in 1:n){
  sliced <- df_estimates[i,]
  W_sliced <- W_df[i+start - 1, 2]
  
  rect(xleft = (sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_95, ybottom = sliced$upper_95,
       col = c_mid, border = NA
  )
  
  rect(xleft =(sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_50, ybottom = sliced$upper_50,
       col = c_mid_highlight, border = NA
  )
  
  points(x=sliced$X, y=W_sliced, pch=1)
}

