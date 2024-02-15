library(tidyverse) # organization tools
library(lubridate) # prep survey date data 
library(bipartite) # calculate species interaction metrics

process_raw_data <- function(min_unique_detections) {

  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
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
    filter(!SPECIES %in% c("Apis mellifera", "undetermined", "undetermined/unconfirmed ID"))  %>%
    
    # Filter by SPECIES
    # Remove some others until resolved for final analyses
    filter(!SPECIES %in% c("Bombus sp.","Eupeodes sp.", "Sphaerophoria sp."))  %>%
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
    # for now we will filter out survey round 7 in year 2 (6 / 18 sites visited a 7th time)
    filter(SAMPLING_ROUND <7)
  
  # number of unique species
  (length(unique(mydata_filtered$SPECIES)))
  
  # easy table summary of total detections per species
  detections_species <- mydata_filtered %>%
    group_by(SPECIES) %>%
    count() %>%
    rename(total_detections = n) 
  
  # remove species detected fewer than min unique occassions (site/year/visit) for figure
  mydata_filtered_min_unique <- mydata_filtered %>%  
    group_by(SPECIES, SITE, YEAR, SAMPLING_ROUND) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(SPECIES) %>%
    add_tally() %>%
    rename(unique_detections = n) %>%
    filter(unique_detections >= min_unique_detections)
  
  species_min_or_more <- unique(mydata_filtered_min_unique$SPECIES)
  temp <- mydata_filtered %>%
    filter(SPECIES %in% species_min_or_more)
  
  # total detections all species
  ggplot(mydata_filtered, aes(x=fct_infreq(SPECIES))) +
    geom_bar(stat = "count") +
    labs(x = "", y = "Total detections") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 11)) +
    ggtitle("Total detections per species")
  
  # total detections (species detected 2 or more unique occasions)
  ggplot(temp, aes(x=fct_infreq(SPECIES))) +
    geom_bar(stat = "count") +
    labs(x = "", y = "Total detections") +
    theme(axis.text.x = element_text(size = 11, angle = 65, vjust = 1, hjust=1.1),
          axis.text.y = element_text(size = 11)) +
    ggtitle(paste0("Total detections per species (for species detected ", min_unique_detections, " or more unique visits)"))
  
  # unique detections (species detected 2 or more unique occasions)
  ggplot(mydata_filtered_min_unique, aes(x=fct_infreq(SPECIES))) +
    geom_bar(stat = "count") +
    labs(x = "", y = "Unique detections") +
    theme(axis.text.x = element_text(size = 11, angle = 65, vjust = 1, hjust=1.1),
          axis.text.y = element_text(size = 11)) +
    ggtitle(paste0("Unique site/year/visit detections per species (for species detected ", min_unique_detections, " or more unique visits)"))
  
  
  ## --------------------------------------------------
  ## Transform into detection non detection data
  
  df <- mydata_filtered %>%
    
    # convert to binary detections
    group_by(SPECIES, SITE, YEAR, SAMPLING_ROUND) %>%
    slice(1) %>%
    ungroup() %>%
    
    # remove species detected on < min_unique_detections unique occassions
    group_by(SPECIES) %>%
    add_tally() %>%
    rename(unique_detections = n) %>%
    filter(unique_detections >= min_unique_detections) %>%
    ungroup() %>%
    
    # change SITE to factor for left_join with HABITAT_CATEGORY created below
    mutate(SITE = as.factor(SITE))
  
  # add habitat_category of sites
  # make a new dataframe with site names and cats
  SITE <- levels(as.factor(df$SITE))
  HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                        0,1,1,1,0,0,0,1,0) # and 9 meadow sites
  
  habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY))
  
  # join HABITAT_CATEGORY of sites with df by SITE
  df_joined <- left_join(df, habitats_df, by = "SITE")
  
  ## Get unique species and sites
  # create an alphabetized list of all species encountered across all sites*intervals*visits
  species_list <- df_joined %>%
    group_by(SPECIES) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(SPECIES) # extract species names column as vector
  
  # create an alphabetized list of all sites
  site_list <- df_joined %>%
    group_by(SITE) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(SITE) # extract species names column as vector
  
  # get vectors of species, sites, intervals, and visits 
  species_vector <- species_list %>%
    pull(SPECIES)
  
  site_vector <- site_list %>%
    pull(SITE)
  
  year_vector <- as.vector(levels(as.factor(df$YEAR)))
  
  visit_vector <- as.vector(levels(as.factor(df$SAMPLING_ROUND)))
  
  # find study dimensions from the pooled data
  n_species <- length(species_vector)
  
  n_sites <- length(site_vector)
  
  n_years <- length(levels(as.factor(df$YEAR)))
  
  n_visits <- length(levels(as.factor(df$SAMPLING_ROUND)))
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V
  
  df_joined <- df_joined %>%
    mutate(y = 1) %>%
    arrange(SAMPLING_ROUND, YEAR, SITE, SPECIES) %>% # arrange opposite to array fill 
    # need numeric versions of species, site, year, and visit
    mutate(SPECIES_numeric = as.numeric(as.factor(SPECIES)),
           SITE_numeric = as.numeric(as.factor(SITE)),
           YEAR_numeric = as.numeric(as.factor(YEAR)),
           SAMPLING_ROUND_numeric = as.numeric(as.factor(SAMPLING_ROUND)))
  
  # prepare array and prefill with NA's
  V <- array(data = 0, dim = c(n_species, n_sites, n_years, n_visits))
  
  # fill array with detection / nondetection data
  for(i in 1:nrow(df_joined)){
    V[df_joined$SPECIES_numeric[i], df_joined$SITE_numeric[i], 
      df_joined$YEAR_numeric[i], df_joined$SAMPLING_ROUND_numeric[i]] <- df_joined$y[i]
  }
  
  array_fill_check <- sum(V) - nrow(df_joined) # quick check that they're identical # should be 0
  
  if(array_fill_check != 0){
    print("WARNING: filled array is inconsistent with number of detections in the 2D datasheet!")
  }
  
  ## --------------------------------------------------
  ## make an matrix (site by year) containing the number of total surveys that occurred
  # this will allow us to have a variable likelihood function that takes 7 visits for site X year combos visited 7 times
  # or 6 visits for site X year combos visited only 6 times (most of the combos are 6 annual visits)
  
  temp <- mydata %>%
    group_by(SITE, YEAR, SAMPLING_ROUND) %>%
    slice(1) %>%
    ungroup() %>%
    select(SITE, YEAR, SAMPLING_ROUND) %>%
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND))))  %>%
    group_by(SITE, YEAR) %>%
    mutate(num_surveys = max(SAMPLING_ROUND)) %>%
    slice(1) %>%
    ungroup() %>%
    
    pivot_wider(., names_from = YEAR, values_from = num_surveys) 
    
  site_year_visit_count <- as.matrix(temp[1:n_sites, 3:(3 + n_years - 1)])
  
  rm(temp)
  
  ## --------------------------------------------------
  ## make an array of the scaled visit dates which will be used as a detection covariate
  
  # need to get all dates from all sites not just those with detections
  # so we need to start from the beginning to pull this information back out
  visits <- mydata %>%
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
    group_by(SAMPLING_ROUND, YEAR, SITE) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    select(SITE, YEAR, SAMPLING_ROUND, DATE) %>%
    
    # need to convert dates to julian day of year
    mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
    # need to centre and scale the dates as a predictor variable (z-score)
    mutate(SCALED_DATE = center_scale(DATE))
  
  # initialize a visit array of dimensions n_sites*n_years*n_visits
  scaled_date_array <- array(NA, dim =c(n_sites, n_years, n_visits))
  
  for(i in 1:n_visits){
    
    temp <- visits %>%
      filter(SAMPLING_ROUND == i) %>%
      select(SITE, YEAR, SCALED_DATE) %>%
      pivot_wider(names_from = YEAR, values_from = SCALED_DATE)
    
    if(i %in% 1:6){ # if there were six visits in the year:
      scaled_date_array[,,i] <- as.matrix(temp[,2:(n_years + 1)])
    
      } else{ # we didnt survey all sites 7 times in all years so we need to handle differently
      
      temp2 <- left_join(site_list, temp) 
      
      # stan can't handle NA's
      # instead we will pick an unfeasible number so that we don't get confused
      # and then later tell stan to just pass over these values (do not evaluate likelihood)
      temp2[is.na(temp2)] <- -99
      
      temp3 <- pull(temp2[,2])
      
      scaled_date_array[,2,i] <- temp3
      scaled_date_array[,1,i] <- rep(-99, times=n_sites)
      scaled_date_array[,3,i] <- rep(-99, times=n_sites)
    }
    
  }
  
  if(length(which(is.na(scaled_date_array))) != 0){
    print("WARNING: missing survey dates for 1 or more site/year/visit combinations!")
  }
  
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
  ## Get site X year flowering plant abundance data (herbaceous plants in quadrats)
  
  # read data
  plant_data <- read.csv("./data/flower_resources_quadrats.csv")
  
  # perform some initial filters on the unfinished prelim data
  plant_data <- plant_data %>% 
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
    
    # NA's indicate NO FLOWERS WERE IN BLOOM 
    mutate(NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0))
  
  # read pollinator data
  plants_visited <- mydata_filtered %>%
    group_by(PLANT_NETTED_FROM_SCI_NAME) %>%
    add_tally() %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(PLANT_NETTED_FROM_SCI_NAME, n) %>% 
    mutate(log_n = log(n)) %>%
    filter(PLANT_NETTED_FROM_SCI_NAME != "") %>%
    mutate(PLANT_NETTED_FROM_SCI_NAME = fct_reorder(PLANT_NETTED_FROM_SCI_NAME, desc(n))) %>%
    filter(n >= 5)
  
  ggplot(plants_visited, aes(x=PLANT_NETTED_FROM_SCI_NAME, y=log_n)) +
    geom_col() +
    labs(x = "Plant species", y="log(number of detected interactions)") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    ggtitle("Total interactions per species")
  
  # create an alphabetized list of all species encountered across all sites*intervals*visits
  plant_species_list_reduced <- plants_visited %>%
    select(PLANT_NETTED_FROM_SCI_NAME) # extract species names column as vector
  
  # now filter out plants that were never visited by pollinators (considered not to be of high value for pollinators)
  plant_data_subset <- plant_data %>%
    filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME)
  
  plant_abundance_df <- plant_data_subset %>%
    
    # calculate abundance per site visit
    group_by(SITE, YEAR, SAMPLING_ROUND) %>%
    mutate(total_flower_abundance = sum(NUM_FLORAL_UNITS),
           log_total_flower_abundance = log(total_flower_abundance + 1)) %>%
    
    # take one row per site visit
    slice(1) %>%
    ungroup() %>%
    
    # there could be multiple ways to group and calculate plant data
    # here grouping by average abundance per year
    group_by(SITE, YEAR) %>%
    mutate(mean_annual_plant_abundance = mean(log_total_flower_abundance)) %>%
    slice(1) %>%
    ungroup() %>%
  
    # scale the variable
    mutate(plant_abundance_scaled = center_scale(mean_annual_plant_abundance))
    
  original_herb_abundance <- plant_abundance_df$mean_annual_plant_abundance
    
  plant_abundance_df <- plant_abundance_df %>%
    select(SITE, YEAR, plant_abundance_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = plant_abundance_scaled)
  
  herbaceous_flowers_scaled <- as.matrix(plant_abundance_df[2:(n_years+1)])
  
  ## --------------------------------------------------
  ## Get site X year flowering plant abundance data (woody plants within survey area)
  
  # read data
  woody_plant_data <- read.csv("./data/flower_resources_woody.csv")
  
  # perform some initial filters on the unfinished prelim data
  woody_plant_data <- woody_plant_data %>% 
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
    
    # NA's indicate NO FLOWERS WERE IN BLOOM 
    mutate(SPECIES = replace_na(SPECIES, "No woody flowers"),
           NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0)) %>%

    filter(SHRUB_OR_TREE == "y") 
  
  site_visits <- as.data.frame(cbind(
    rep(site_vector, each=3*6), rep(1:6, times=18*3), rep(1:3, each = 6, times=3), 
    rep("no flowering species", 324), 
    rep(0, 324))) %>%
    rename("SITE" = "V1",
           "SAMPLING_ROUND" = "V2",
           "YEAR" = "V3",
           "SPECIES" = "V4",
           "NUM_FLORAL_UNITS" = "V5") %>%
    mutate(YEAR = as.integer(YEAR),
           SAMPLING_ROUND = as.integer(SAMPLING_ROUND), 
           NUM_FLORAL_UNITS = as.integer(NUM_FLORAL_UNITS))
  
  # now filter out plants that were never visited by pollinators (considered not to be of high value for pollinators)
  # we already created this list above when looking at the quadrat data
  woody_plant_data_subset <- woody_plant_data %>%
    filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME)
  
  # join all possible site visits back in to add zeros
  woody_plant_data_subset <- full_join(woody_plant_data_subset, site_visits) 
  
  # get mean per site X year
  woody_plant_data_subset <- woody_plant_data_subset %>%
    mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%
    group_by(SITE, YEAR) %>%
    mutate(mean_annual_woody_plant_abundance = mean(log_NUM_FLORAL_UNITS)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(woody_plant_abundance_scaled = center_scale(mean_annual_woody_plant_abundance)) 
  
  original_woody_abundance <- woody_plant_data_subset$mean_annual_woody_plant_abundance
  
  woody_plant_data_subset <- woody_plant_data_subset %>%
    select(SITE, YEAR, woody_plant_abundance_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = woody_plant_abundance_scaled)
  
  woody_flowers_scaled <- as.matrix(woody_plant_data_subset[2:(n_years+1)])
  
  ## --------------------------------------------------
  ## Get visit specific flowering plant abundance data (for detection covariate)
  # survey_flower_abundance should be an array of size [1:n_sites, 1:n_years, 1:n_visits]
  
  woody_plant_data_subset <- woody_plant_data %>%
    filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME)
  
  # join all possible site visits back in to add zeros
  woody_plant_data_subset <- full_join(woody_plant_data_subset, site_visits)
  
  # get site X year X sampling round
  woody_plant_data_subset <- woody_plant_data_subset %>%
    mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%
    group_by(SITE, YEAR, SAMPLING_ROUND) %>%
    mutate(survey_woody_plant_abundance = sum(log_NUM_FLORAL_UNITS)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(survey_woody_plant_abundance_scaled = center_scale(survey_woody_plant_abundance)) %>%
    select(SITE, YEAR, SAMPLING_ROUND, survey_woody_plant_abundance_scaled) %>%
    mutate(SITE_NUMBER = rep(1:n_sites, each = n_years*n_visits))
    
  woody_flowers_by_survey = array(NA, dim = c(n_sites,n_years,n_visits))
  
  for(i in 1:nrow(woody_plant_data_subset)){

    woody_flowers_by_survey[woody_plant_data_subset$SITE_NUMBER[i], 
                          woody_plant_data_subset$YEAR[i],
                          woody_plant_data_subset$SAMPLING_ROUND[i]] <- 
      woody_plant_data_subset$survey_woody_plant_abundance_scaled[i]

  }
  
  # now do the same for the herbaceous plants
  herabceous_plant_data_subset <- full_join(plant_data_subset, site_visits)

  # get site X year X sampling round
  herabceous_plant_data_subset <- herabceous_plant_data_subset %>%
    mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%
    group_by(SITE, YEAR, SAMPLING_ROUND) %>%
    mutate(survey_herbaceous_plant_abundance = sum(log_NUM_FLORAL_UNITS)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(survey_herbaceous_plant_abundance_scaled = center_scale(survey_herbaceous_plant_abundance)) %>%
    select(SITE, YEAR, SAMPLING_ROUND, survey_herbaceous_plant_abundance_scaled) %>%
    mutate(SITE_NUMBER = rep(1:n_sites, each = n_years*n_visits))
  
  herbaceous_flowers_by_survey = array(NA, dim = c(n_sites,n_years,n_visits))
  
  for(i in 1:nrow(herabceous_plant_data_subset)){
    
    herbaceous_flowers_by_survey[herabceous_plant_data_subset$SITE_NUMBER[i], 
                          herabceous_plant_data_subset$YEAR[i],
                          herabceous_plant_data_subset$SAMPLING_ROUND[i]] <- 
      herabceous_plant_data_subset$survey_herbaceous_plant_abundance_scaled[i]
    
  }
  
  # now make a composite average deviation from average abundance
  # across woody and herbaceous plants (in case we want to only model single parameter for effects of all flowers)
  flowers_any_by_survey = array(NA, dim = c(n_sites,n_years,n_visits))
  
  flowers_any_by_survey = (woody_flowers_by_survey + herbaceous_flowers_by_survey) / 2
  
  ## --------------------------------------------------
  ## Return stuff
  return(list(
    
    n_species = n_species, ## number of species
    n_sites = n_sites, ## number of sites
    n_years = n_years, ## number of occupancy intervals (years)
    n_visits = n_visits, ## number of samples per year
    V = V, # detection array
    species = species_vector, # alphabetized species names
    sites = site_vector, # alphabetized site names
    years = year_vector, # ordered vector of intervals
    visits = visit_vector, # ordered vector of visits
    site_year_visit_count = site_year_visit_count, # matrix containing number of visits per site X year 
    date_scaled = scaled_date_array, # detection covariate (array of scaled julian date of visit)
    habitat_category = HABITAT_CATEGORY, # occupancy and detection covariate (sites mowed or meadow)
    species_interaction_metrics = species_interaction_metrics,
    herbaceous_flowers_scaled = herbaceous_flowers_scaled,
    woody_flowers_scaled = woody_flowers_scaled,
    herbaceous_flowers_by_survey = herbaceous_flowers_by_survey,
    woody_flowers_by_survey = woody_flowers_by_survey,
    flowers_any_by_survey = flowers_any_by_survey,
    
    original_herb_abundance = original_herb_abundance,
    original_woody_abundance = original_woody_abundance
    
  ))
  
}
