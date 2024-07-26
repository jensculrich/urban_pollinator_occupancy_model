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
  # remove ~2 dozen rows where we accidentally printed an extra label 
  # (i.e., field assistant listed more species interacting with a plant at a site/visit than were collected)
  filter(!SPECIES %in% c("Apis mellifera", "extra_label"))  %>%
  
  # Filter by SPECIES
  # Remove some others until resolved for final analyses
  # Two un'id bombus were failed to be collected and brought back to lab for ID
  # one Eupeodes had an unconfident ID (looks intermediate between two species) so we are excluding it
  # No male Osmia were given an ID and were left simply as Osmia sp.
  filter(!SPECIES %in% c("Bombus sp.","Eupeodes sp.", "Osmia sp."))  %>%
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
  # for now we will filter out survey round 7 in year 2 (6 / 18 sites visited a 7th time)
  filter(SAMPLING_ROUND < 7)

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
## Species interaction frequency data to be used for model

# group by species, and calculate number of plants per species
interactions_df <- mydata_filtered %>% 
  rename("pollinator_species" = "SPECIES",
         "plant_species" = "PLANT_NETTED_FROM_SCI_NAME")  %>% 
  select(pollinator_species, plant_species) %>%
  group_by(pollinator_species, plant_species) %>% 
  add_tally() %>%
  rename("connection_strength" = "n") %>%
  slice(1) %>%
  ungroup() 

# get all of the plant species that pollinators were recorded to interact with
plant_species <- mydata_filtered %>%
  group_by(PLANT_NETTED_FROM_SCI_NAME) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(PLANT_NETTED_FROM_SCI_NAME) # extract species names column as vector

# plant names
plant_species_vector <- plant_species %>%
  pull(PLANT_NETTED_FROM_SCI_NAME)

# number of plant names
n_plant_species <- nrow(plant_species)

# get all possible interactions so that we can fill undetected inteactions with zeroes
all_possible_interactions <- as.data.frame(cbind(rep(species_vector, each = n_plant_species),
                                                 rep(plant_species_vector, times = n_species))) %>%
  rename("pollinator_species" = "V1",
         "plant_species" = "V2")
# join all possible interactions with observed interactions
interactions_df <- left_join(all_possible_interactions, interactions_df) 
# replace missing interactions with zeroes instead of NA
interactions_df$connection_strength <- interactions_df$connection_strength %>% replace_na(0)

# prep an interaction frequency matrix
A <- matrix(data = 0, nrow = n_species, ncol=n_plant_species,
            dimnames = list(species_vector, plant_species_vector))

# fill matrix with observed interaction strengths
for(i in 1:nrow(interactions_df)){
  A[interactions_df$pollinator_species[i], interactions_df$plant_species[i]] <- 
    interactions_df$connection_strength[i]
}

# calculate species level interaction metrics
species_interaction_metrics <- bipartite::specieslevel(A, level="lower", index=c("degree", "normalised degree", "d"), PDI.normalise=F)

species_interaction_metrics <- species_interaction_metrics %>%
  mutate(degree_scaled = center_scale(degree),
         normalised.degree_scaled = center_scale(normalised.degree),
         d_scaled = center_scale(d))

## --------------------------------------------------
## Species interactions SUPPLEMENTED WITH ELIZABETH ELLE'S INTERACTION DATA

# group by species, and calculate number of plants per species
interactions_df2 <- mydata_filtered %>% 
  rename("pollinator_species" = "SPECIES",
         "plant_species" = "PLANT_NETTED_FROM_SCI_NAME")  %>% 
  select(pollinator_species, plant_species)

extra_interactions <- read.csv("./data/elle_data_2/reduced_pollinator_plant.csv") %>%
  mutate(pollinator_species = paste(pollinator_genus, pollinator_species, sep = " "))

# prep the data so it is formatted similarly

# one thing to note is that we can't have correspondance between morphospecies
# (because there's no way to match up who's who across the two datasets)
# I use the term "morph 1" while Elle uses the term "sp 1", so they should never match
# but unfortunately this means we can't bolster our interaction detections for these pollinators

# rename some species used in the Elle dataset
# either where they used a different naming concept
# or where they used a narrower taxonomic scope than I did
extra_interactions <- extra_interactions %>%
  mutate(pollinator_species = gsub('Bombus bifarius', 'Bombus vancouverensis', pollinator_species),
         pollinator_species = gsub('Bombus californicus', 'Bombus fervidus (californicus)', pollinator_species),
         pollinator_species = gsub('Coelioxys porterae', 'Coelioxys sp.', pollinator_species),
         pollinator_species = gsub('Eumerus sp', 'Eumerus sp.', pollinator_species),
         pollinator_species = gsub('Eumerus narcissi', 'Eumerus sp.', pollinator_species),
         pollinator_species = gsub('Eupeodes lapponicus', 'Lapposyrphus lapponicus', pollinator_species),
         pollinator_species = gsub('Lasioglossum titusi', 'Lasioglossum titusii', pollinator_species),
         pollinator_species = gsub('Melissodes microsticta', 'Melissodes microstictus', pollinator_species),
         pollinator_species = gsub('Allograpta micrura', 'Fazia macrura', pollinator_species),
         pollinator_species = gsub('Merodon Merodon equestris', 'Merodon equestris', pollinator_species),
  )


# rename some plant species used in the Elle dataset
# either where they used a different naming concept
# or where they used a narrower taxonomic scope than I did
extra_interactions <- extra_interactions %>%
  mutate(plant_species = gsub('^Achillea$', 'Achillea_millefolium', plant_species),
         plant_species = gsub('Buddleja davidii', 'Buddleja', plant_species),
         plant_species = gsub('Buddleja globosa', 'Buddleja', plant_species),
         plant_species = gsub('Buddleja spp.', 'Buddleja', plant_species),
         plant_species = gsub('Brassica - yellow', 'Brassica_sp.', plant_species),
         plant_species = gsub('^Centaurea (garden)$', 'Centaurea', plant_species),
         plant_species = gsub('Cerastium (garden)', 'Cerastium arvense', plant_species),
         plant_species = gsub('^Centaurea sp.$', 'Centaurea', plant_species),
         plant_species = gsub('Chamaenerion angustifolium', 'Chamerion angustifolium', plant_species),
         plant_species = gsub('^Cirsium vulgare$', 'Cirsium vulgaris', plant_species),
         plant_species = gsub('Clematis', 'Clematis vitalba', plant_species),
         plant_species = gsub('^Coreopsis$', 'Coreopsis_lanceolata', plant_species),
         plant_species = gsub('Craetagus douglasii', 'Crataegus douglasii', plant_species),
         plant_species = gsub('^Digitalis$', 'Digitalis purpurea', plant_species),
         plant_species = gsub('^Echinacea$', 'Echinacea purpurea', plant_species),
         plant_species = gsub('^Fragaria chiloensis$', 'Fragaria', plant_species),
         plant_species = gsub('^Fragaria virginiana$', 'Fragaria', plant_species),
         plant_species = gsub('^Hypericum$', 'Hypericum_perforatum', plant_species),
         plant_species = gsub('^Iris$', 'Iris sp.', plant_species),
         plant_species = gsub('^Lamium$', 'Lamium purpureum', plant_species),
         plant_species = gsub('^Lamium purpureum (invasive)$', 'Lamium purpureum', plant_species),
         plant_species = gsub('^Limnanthes$', 'Limnanthes douglasii', plant_species),
         plant_species = gsub('^Limnanthes douglassi$', 'Limnanthes douglasii', plant_species),
         plant_species = gsub('Malus fusca', 'Malus', plant_species),
         plant_species = gsub('Papaver sp', 'Papaver_sp.', plant_species),
         plant_species = gsub('Rubus discolor', 'Rubus_armeniacus', plant_species),
         plant_species = gsub('Rubus parviflora', 'Rubus_parviflorus', plant_species),
         plant_species = gsub('Sedum sp. spathulifolium', 'sp.', plant_species),
         plant_species = gsub('Sedum spurium', 'sp.', plant_species),
         plant_species = gsub('Trifolium pratens', 'Trifolium_pratense', plant_species),
         plant_species = gsub('Sedum', 'Sedum_sp.', plant_species),
         plant_species = gsub('Sorbus aucuparia', 'Sorbus', plant_species),
         plant_species = gsub('Sorbus hybrida', 'Sorbus', plant_species),
         plant_species = gsub('^Spiraea sp. douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spirea sp. douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spirea douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spiraea$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Rhododendron$', 'Rhododendron_macrophyllum', plant_species),
         plant_species = gsub('^Rudbeckia$', 'Rudbeckia_hirta', plant_species),
         plant_species = gsub('Rudbeckia fulgida', 'Rudbeckia_hirta', plant_species)
  )


# and reduce the resolution on some of my obs
interactions_df2 <- interactions_df2 %>% 
  mutate(plant_species = gsub('Buddleja davidii', 'Buddleja', plant_species),
         plant_species = gsub('Campanula rapunculoides', 'Campanula', plant_species),
         plant_species = gsub('^Centaurea cyanus$', 'Centaurea', plant_species),
         plant_species = gsub('Dianthus barbatus', 'Dianthus', plant_species),
         plant_species = gsub('Hyacinthoides non-scripta', 'Hyacinthoides', plant_species),
         plant_species = gsub('Malus domestica', 'Malus', plant_species),
         plant_species = gsub('Malus fusca', 'Malus', plant_species),
         plant_species = gsub('Salvia nemorosa', 'Salvia', plant_species),
         plant_species = gsub('Salvia sp.', 'Salvia', plant_species),
         plant_species = gsub('Sorbus aucuparia', 'Sorbus', plant_species),
         plant_species = gsub('Sorbus hybrida', 'Sorbus', plant_species),
         plant_species = gsub('Spiraea douglasii', 'Spiraea sp.', plant_species)
  )

# group by species, and calculate number of plants per species
extra_interactions_df <- extra_interactions %>%
  #filter(ecosystem != "Shrub-Steppe") %>% # interactions not comparable from different ecosystem?
  # only use th interactions from BC (not Oregon and Washington)
  filter(latitude > 48) %>%
  # replace "_" with " " to match my data
  mutate(pollinator_species = gsub('_', ' ', pollinator_species),
         plant_species = gsub('_', ' ', plant_species)) %>%
  select(pollinator_species, plant_species) %>%
  filter(pollinator_species %in% species_vector) %>%
  rbind(., interactions_df2) %>%
  group_by(pollinator_species, plant_species) %>% 
  add_tally() %>%
  rename("connection_strength" = "n") %>%
  slice(1) %>%
  ungroup() 

# get number of total interactions per plant species
interactions_per_plant_species <- extra_interactions_df %>%
  group_by(plant_species) %>%
  mutate(total_interactions = sum(connection_strength)) %>%
  slice(1) %>%
  ungroup %>%
  select(plant_species, total_interactions)

# get all of the plant species that pollinators were recorded to interact with
plant_species_orig <- mydata_filtered %>%
  group_by(PLANT_NETTED_FROM_SCI_NAME) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(PLANT_NETTED_FROM_SCI_NAME) %>% # extract species names column as vector
  rename("plant_species" = "PLANT_NETTED_FROM_SCI_NAME")

# get all of the plant species that pollinators were recorded to interact with
plant_species_elle <- extra_interactions_df %>%
  group_by(plant_species) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(plant_species) # extract species names column as vector

# plant names
plant_species_vector_w_extra_data <- plant_species_orig %>%
  rbind(., plant_species_elle) %>%
  group_by(plant_species) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  pull(plant_species)

# number of plant names
n_plant_species_w_extra_data <- length(plant_species_vector_w_extra_data)

# get all possible interactions so that we can fill undetected inteactions with zeroes
all_possible_interactions2 <- as.data.frame(cbind(rep(species_vector, each = n_plant_species_w_extra_data),
                                                  rep(plant_species_vector_w_extra_data, times = n_species))) %>%
  rename("pollinator_species" = "V1",
         "plant_species" = "V2")
# join all possible interactions with observed interactions
interactions_df2 <- left_join(all_possible_interactions2, extra_interactions_df) 
# replace missing interactions with zeroes instead of NA
interactions_df2$connection_strength <- interactions_df2$connection_strength %>% replace_na(0)

# prep an interaction frequency matrix
A2 <- matrix(data = 0, nrow = n_species, ncol=n_plant_species_w_extra_data,
             dimnames = list(species_vector, plant_species_vector_w_extra_data))

# fill matrix with observed interaction strengths
for(i in 1:nrow(interactions_df2)){
  A2[interactions_df2$pollinator_species[i], interactions_df2$plant_species[i]] <- 
    interactions_df2$connection_strength[i]
}

# calculate species level interaction metrics
species_interaction_metrics2 <- bipartite::specieslevel(A2, level="lower", index=c("degree", "normalised degree", "d"), PDI.normalise=F)

species_interaction_metrics2 <- species_interaction_metrics2 %>%
  mutate(degree_scaled = center_scale(degree),
         normalised.degree_scaled = center_scale(normalised.degree),
         d_scaled = center_scale(d)) %>%
  rename("degree_supplemented" = "degree",
         "normalised.degree_supplemented" = "normalised.degree",
         "d_supplemented" = "d",
         "degree_scaled_supplemented" = "degree_scaled",
         "normalised.degree_scaled_supplemented" = "normalised.degree_scaled",
         "d_scaled_supplemented" = "d_scaled")

species_interaction_metrics2 <- cbind(species_interaction_metrics, species_interaction_metrics2)

species_interaction_metrics <- species_interaction_metrics2

cor <- cor(species_interaction_metrics$d_scaled, species_interaction_metrics$d_scaled_supplemented)
par(mar = c(5, 5, 2, 2))
plot(species_interaction_metrics$d_scaled, species_interaction_metrics$d_scaled_supplemented,
     xlab="specialization (my data alone)", ylab="specialization (with Elle lab interactions)",)
text(x = 3, y = -1, # Coordinates
     label = paste0("Pearson's R = ", signif(cor, 3)))


## --------------------------------------------------
## Species interactions at the genus level SUPPLEMENTED WITH ELIZABETH ELLE'S INTERACTION DATA

# group by species, and calculate number of plants per species
interactions_df2 <- mydata_filtered %>% 
  rename("pollinator_species" = "SPECIES",
         "plant_species" = "PLANT_NETTED_FROM_SCI_NAME")  %>% 
  select(pollinator_species, plant_species)

extra_interactions <- read.csv("./data/elle_data_2/reduced_pollinator_plant.csv") %>%
  mutate(pollinator_species = paste(pollinator_genus, pollinator_species, sep = " "))

# prep the data so it is formatted similarly

# one thing to note is that we can't have correspondance between morphospecies
# (because there's no way to match up who's who across the two datasets)
# I use the term "morph 1" while Elle uses the term "sp 1", so they should never match
# but unfortunately this means we can't bolster our interaction detections for these pollinators

# rename some species used in the Elle dataset
# either where they used a different naming concept
# or where they used a narrower taxonomic scope than I did
extra_interactions <- extra_interactions %>%
  mutate(pollinator_species = gsub('Bombus bifarius', 'Bombus vancouverensis', pollinator_species),
         pollinator_species = gsub('Bombus californicus', 'Bombus fervidus (californicus)', pollinator_species),
         pollinator_species = gsub('Coelioxys porterae', 'Coelioxys sp.', pollinator_species),
         pollinator_species = gsub('Eumerus sp', 'Eumerus sp.', pollinator_species),
         pollinator_species = gsub('Eumerus narcissi', 'Eumerus sp.', pollinator_species),
         pollinator_species = gsub('Eupeodes lapponicus', 'Lapposyrphus lapponicus', pollinator_species),
         pollinator_species = gsub('Lasioglossum titusi', 'Lasioglossum titusii', pollinator_species),
         pollinator_species = gsub('Melissodes microsticta', 'Melissodes microstictus', pollinator_species),
         pollinator_species = gsub('Allograpta micrura', 'Fazia macrura', pollinator_species),
         pollinator_species = gsub('Merodon Merodon equestris', 'Merodon equestris', pollinator_species),
  )


# rename some plant species used in the Elle dataset
# either where they used a different naming concept
# or where they used a narrower taxonomic scope than I did
extra_interactions <- extra_interactions %>%
  mutate(plant_species = gsub('^Achillea$', 'Achillea_millefolium', plant_species),
         plant_species = gsub('Buddleja davidii', 'Buddleja', plant_species),
         plant_species = gsub('Buddleja globosa', 'Buddleja', plant_species),
         plant_species = gsub('Buddleja spp.', 'Buddleja', plant_species),
         plant_species = gsub('Brassica - yellow', 'Brassica_sp.', plant_species),
         plant_species = gsub('^Centaurea (garden)$', 'Centaurea', plant_species),
         plant_species = gsub('Cerastium (garden)', 'Cerastium arvense', plant_species),
         plant_species = gsub('^Centaurea sp.$', 'Centaurea', plant_species),
         plant_species = gsub('Chamaenerion angustifolium', 'Chamerion angustifolium', plant_species),
         plant_species = gsub('^Cirsium vulgare$', 'Cirsium vulgaris', plant_species),
         plant_species = gsub('Clematis', 'Clematis vitalba', plant_species),
         plant_species = gsub('^Coreopsis$', 'Coreopsis_lanceolata', plant_species),
         plant_species = gsub('Craetagus douglasii', 'Crataegus douglasii', plant_species),
         plant_species = gsub('^Digitalis$', 'Digitalis purpurea', plant_species),
         plant_species = gsub('^Echinacea$', 'Echinacea purpurea', plant_species),
         plant_species = gsub('^Fragaria chiloensis$', 'Fragaria', plant_species),
         plant_species = gsub('^Fragaria virginiana$', 'Fragaria', plant_species),
         plant_species = gsub('^Hypericum$', 'Hypericum_perforatum', plant_species),
         plant_species = gsub('^Iris$', 'Iris sp.', plant_species),
         plant_species = gsub('^Lamium$', 'Lamium purpureum', plant_species),
         plant_species = gsub('^Lamium purpureum (invasive)$', 'Lamium purpureum', plant_species),
         plant_species = gsub('^Limnanthes$', 'Limnanthes douglasii', plant_species),
         plant_species = gsub('^Limnanthes douglassi$', 'Limnanthes douglasii', plant_species),
         plant_species = gsub('Malus fusca', 'Malus', plant_species),
         plant_species = gsub('Papaver sp', 'Papaver_sp.', plant_species),
         plant_species = gsub('Rubus discolor', 'Rubus_armeniacus', plant_species),
         plant_species = gsub('Rubus parviflora', 'Rubus_parviflorus', plant_species),
         plant_species = gsub('Sedum sp. spathulifolium', 'sp.', plant_species),
         plant_species = gsub('Sedum spurium', 'sp.', plant_species),
         plant_species = gsub('Trifolium pratens', 'Trifolium_pratense', plant_species),
         plant_species = gsub('Sedum', 'Sedum_sp.', plant_species),
         plant_species = gsub('Sorbus aucuparia', 'Sorbus', plant_species),
         plant_species = gsub('Sorbus hybrida', 'Sorbus', plant_species),
         plant_species = gsub('^Spiraea sp. douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spirea sp. douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spirea douglasii ssp. Douglasii$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Spiraea$', 'Spiraea_sp.', plant_species),
         plant_species = gsub('^Rhododendron$', 'Rhododendron_macrophyllum', plant_species),
         plant_species = gsub('^Rudbeckia$', 'Rudbeckia_hirta', plant_species),
         plant_species = gsub('Rudbeckia fulgida', 'Rudbeckia_hirta', plant_species)
  )


# and reduce the resolution on some of my obs
interactions_df2 <- interactions_df2 %>% 
  mutate(plant_species = gsub('Buddleja davidii', 'Buddleja', plant_species),
         plant_species = gsub('Campanula rapunculoides', 'Campanula', plant_species),
         plant_species = gsub('^Centaurea cyanus$', 'Centaurea', plant_species),
         plant_species = gsub('Dianthus barbatus', 'Dianthus', plant_species),
         plant_species = gsub('Hyacinthoides non-scripta', 'Hyacinthoides', plant_species),
         plant_species = gsub('Malus domestica', 'Malus', plant_species),
         plant_species = gsub('Malus fusca', 'Malus', plant_species),
         plant_species = gsub('Salvia nemorosa', 'Salvia', plant_species),
         plant_species = gsub('Salvia sp.', 'Salvia', plant_species),
         plant_species = gsub('Sorbus aucuparia', 'Sorbus', plant_species),
         plant_species = gsub('Sorbus hybrida', 'Sorbus', plant_species),
         plant_species = gsub('Spiraea douglasii', 'Spiraea sp.', plant_species)
  )

# group by species, and calculate number of plants per species
extra_interactions_df <- extra_interactions %>%
  #filter(ecosystem != "Shrub-Steppe") %>% # interactions not comparable from different ecosystem?
  # only use th interactions from BC (not Oregon and Washington)
  filter(latitude > 48) %>%
  # replace "_" with " " to match my data
  mutate(pollinator_species = gsub('_', ' ', pollinator_species),
         plant_species = gsub('_', ' ', plant_species)) %>%
  select(pollinator_species, plant_species) %>%
  filter(pollinator_species %in% species_vector) %>%
  rbind(., interactions_df2) %>%
  # now add the genus (first word of plant_species)
  mutate(plant_genus = word(plant_species, 1)) %>% 
  group_by(pollinator_species, plant_genus) %>% 
  add_tally() %>%
  rename("connection_strength" = "n") %>%
  slice(1) %>%
  ungroup() %>%
  select(-plant_species)

# get all of the plant genera that pollinators were recorded to interact with
plant_genus_orig <- mydata_filtered %>%
  group_by(PLANT_NETTED_FROM_GENUS) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(PLANT_NETTED_FROM_GENUS) %>% # extract species names column as vector
  rename("plant_genus" = "PLANT_NETTED_FROM_GENUS")

# get all of the plant species that pollinators were recorded to interact with
plant_genus_elle <- extra_interactions_df %>%
  group_by(plant_genus) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(plant_genus) # extract species names column as vector

# plant names
plant_genus_vector_w_extra_data <- plant_genus_orig %>%
  rbind(., plant_genus_elle) %>%
  group_by(plant_genus) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  pull(plant_genus)

# number of plant names
n_plant_genus_vector_w_extra_data <- length(plant_genus_vector_w_extra_data)

# get all possible interactions so that we can fill undetected inteactions with zeroes
all_possible_interactions2 <- as.data.frame(cbind(rep(species_vector, each = n_plant_genus_vector_w_extra_data),
                                                  rep(plant_genus_vector_w_extra_data, times = n_species))) %>%
  rename("pollinator_species" = "V1",
         "plant_genus" = "V2")
# join all possible interactions with observed interactions
interactions_df2 <- left_join(all_possible_interactions2, extra_interactions_df) 
# replace missing interactions with zeroes instead of NA
interactions_df2$connection_strength <- interactions_df2$connection_strength %>% replace_na(0)

# prep an interaction frequency matrix
A2 <- matrix(data = 0, nrow = n_species, ncol=n_plant_genus_vector_w_extra_data,
             dimnames = list(species_vector, plant_genus_vector_w_extra_data))

# fill matrix with observed interaction strengths
for(i in 1:nrow(interactions_df2)){
  A2[interactions_df2$pollinator_species[i], interactions_df2$plant_genus[i]] <- 
    interactions_df2$connection_strength[i]
}

# calculate species level interaction metrics
species_interaction_metrics2 <- bipartite::specieslevel(A2, level="lower", index=c("degree", "normalised degree", "d"), PDI.normalise=F)

species_interaction_metrics2 <- species_interaction_metrics2 %>%
  mutate(degree_scaled = center_scale(degree),
         normalised.degree_scaled = center_scale(normalised.degree),
         d_scaled = center_scale(d)) %>%
  rename("degree_supplemented_genus" = "degree",
         "normalised.degree_supplemented_genus" = "normalised.degree",
         "d_supplemented_genus" = "d",
         "degree_scaled_supplemented_genus" = "degree_scaled",
         "normalised.degree_scaled_supplemented_genus" = "normalised.degree_scaled",
         "d_scaled_supplemented_genus" = "d_scaled")

species_interaction_metrics2 <- cbind(species_interaction_metrics, species_interaction_metrics2)

species_interaction_metrics <- species_interaction_metrics2

cor <- cor(species_interaction_metrics$d_scaled, species_interaction_metrics$d_scaled_supplemented_genus)
par(mar = c(5, 5, 2, 2))
plot(species_interaction_metrics$d_scaled, species_interaction_metrics$d_scaled_supplemented_genus,
     xlab="specialization (my data alone)", ylab="specialization at genus level (with Elle lab interactions)",)
text(x = 3, y = -1, # Coordinates
     label = paste0("Pearson's R = ", signif(cor, 3)))
dev.off()


## --------------------------------------------------
## Plot the interactions

plotweb(A2)
visweb(A2)

species_metrics <- specieslevel(A, level="lower", index=c("degree", "normalised degree", "d"), PDI.normalise=F)

hist(species_interaction_metrics$d_supplemented_genus)
hist(species_metrics$degree)

mean(species_interaction_metrics$d_supplemented_genus)
sd(species_interaction_metrics$d_supplemented_genus)

species_interaction_metrics