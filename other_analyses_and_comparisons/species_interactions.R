library(tidyverse)
library(lubridate)
library(bipartite)

## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## --------------------------------------------------
## Read in data 

# read data
mydata <- read.csv("./data/pollinator_data.csv")

# perform some initial filters on the unfinished prelim data
mydata_filtered <- mydata %>% 
  
  # Filter by CLADE. 
  # NA's are specimens that haven't been identified yet. 
  # There are also specimens labelled "Other" if we collected them but they were not bees or hoverflies
  filter(CLADE %in% c("Anthophila", "Syrphidae")) %>%
  
  # Filter by SPECIES
  # Remove honeybees from our preliminary analysis
  filter(!SPECIES %in% c("Apis mellifera","undetermined", "undetermined/unconfirmed ID"))  %>%
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND))))
# %>%
  #filter(!is.na(NO_ID_RESOLUTION))

pollinator_species <- mydata_filtered %>%
  group_by(SPECIES) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SPECIES) # extract species names column as vector

pollinator_species_vector <- pollinator_species %>%
  pull(SPECIES)

plant_species <- mydata_filtered %>%
  group_by(PLANT_NETTED_FROM_SCI_NAME) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(PLANT_NETTED_FROM_SCI_NAME) # extract species names column as vector

plant_species_vector <- plant_species %>%
  pull(PLANT_NETTED_FROM_SCI_NAME)

n_pollinator_species <- nrow(pollinator_species)
n_plant_species <- nrow(plant_species)

## --------------------------------------------------
## Raw number of unique plant species interactions recorded

# group by species, and calculate number of plants per species
unique_plant_interactions <- mydata_filtered %>% 
  group_by(SPECIES, PLANT_NETTED_FROM_SCI_NAME) %>% 
  slice(1) %>%
  ungroup() %>%
  
  group_by(SPECIES) %>%
  add_tally() %>%
  slice(1) %>%
  ungroup() %>%
  
  select(SPECIES, n) %>%
  rename("unique_plant_interactions" = "n")


## --------------------------------------------------
## Visualize network

# group by species, and calculate number of plants per species
test <- mydata_filtered %>% 
  rename("pollinator_species" = "SPECIES",
         "plant_species" = "PLANT_NETTED_FROM_SCI_NAME")  %>% 
  select(pollinator_species, plant_species) %>%
  group_by(pollinator_species, plant_species) %>% 
  add_tally() %>%
  rename("connection_strength" = "n") %>%
  slice(1) %>%
  ungroup() 

all_possible_interactions <- as.data.frame(cbind(rep(pollinator_species_vector, each = n_plant_species),
                                             rep(plant_species_vector, times = n_pollinator_species))) %>%
  rename("pollinator_species" = "V1",
         "plant_species" = "V2")

test2 <- left_join(all_possible_interactions, test) 
  
test2$connection_strength <- test2$connection_strength %>% replace_na(0)

A <- matrix(data = 0, nrow = n_pollinator_species, ncol=n_plant_species,
              dimnames = list(pollinator_species_vector, plant_species_vector))

for(i in 1:nrow(test2)){
  A[test2$pollinator_species[i], test2$plant_species[i]] <- test2$connection_strength[i]
}

plotweb(A)
visweb(A)

species_metrics <- specieslevel(A, level="lower", index=c("degree", "normalised degree", "d"), PDI.normalise=F)

hist(species_metrics$d)
hist(species_metrics$degree)
