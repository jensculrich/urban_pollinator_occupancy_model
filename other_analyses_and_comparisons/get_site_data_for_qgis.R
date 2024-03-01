library(tidyverse)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 1 # >=
my_data <- process_raw_data(min_unique_detections)


## --------------------------------------------------
### Pull out data

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
species_interaction_metrics <- my_data$species_interaction_metrics
d <- species_interaction_metrics$d_scaled
species_names <- my_data$species
site_names <- my_data$sites
herbaceous_flowers_scaled <- my_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_data$woody_flowers_scaled

original_herb_abundance <- my_data$original_herb_abundance
original_woody_abundance <- my_data$original_woody_abundance

# get means across all years by site for summary figure
mean_herb_by_site <- as.data.frame(rowMeans(herbaceous_flowers_scaled))
mean_woody_by_site <- as.data.frame(rowMeans(woody_flowers_scaled))  

# read in site data to join in df
sites <- read.csv("./data/sites.csv", fileEncoding="latin1")

sites <- sites[order(sites$site),]

sites_2 <- cbind(sites, mean_herb_by_site, mean_woody_by_site)

cor.test(sites_2$'rowMeans(herbaceous_flowers_scaled)', sites_2$'rowMeans(woody_flowers_scaled)')

# rename long column names
colnames(sites_2) <- c('site','latitude','longitude', 'category', 'mean_herb_scaled', 'mean_woody_scaled') 

# write the csv file
write.csv(sites_2, "./data/sites_2.csv")
