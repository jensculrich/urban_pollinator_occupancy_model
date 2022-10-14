### Prep real data from urban parks pollinator communities for analysis using occupancy model
# jcu; started oct 13, 2022

library(tidyverse)
library(lubridate)

process_raw_data <- function() {
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  ## --------------------------------------------------
  ## Read data
  df_2021 <- read.csv("./analysis/bee_data_2021.csv")
  df_2022 <- read.csv("./analysis/bee_data_2022.csv")
  
  ## --------------------------------------------------
  ## Prep and clean annual data
  
  # prep 2021 data
  df_2021_prepped <- df_2021 %>%
    filter(CLADE == "Syrphidae") %>% # for now only looking at Syrphids
    group_by(SAMPLING_ROUND, SITE, INSECT_ID) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    mutate(SAMPLING_ROUND_ADJUSTED = SAMPLING_ROUND - 1) %>% # (not using first visits data bc only surveyed for bumblebees)
    mutate(INTERVAL = as.numeric(0)) %>% # year 0
    select(INSECT_ID, SITE, INTERVAL, SAMPLING_ROUND_ADJUSTED, DATE) %>%
    rename("SPECIES" = "INSECT_ID")
  
  # prep 2022 data
  df_2022_prepped <- df_2022 %>%
    filter(CLADE == "Syrphidae") %>% # for now only looking at Syrphids
    group_by(SAMPLING_ROUND, SITE, INSECT_ID) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    mutate(SAMPLING_ROUND_ADJUSTED = SAMPLING_ROUND) %>% 
    mutate(INTERVAL = as.numeric(1)) %>% # year 1
    filter(SAMPLING_ROUND_ADJUSTED != "7") %>% # currently the model can't handle missing data
    # need to find a way to update the model, otherwise we can't use round 7 where only some sites were visited, 
    # and no sites visited in 2021.
    select(INSECT_ID, SITE, INTERVAL, SAMPLING_ROUND_ADJUSTED, DATE) %>%
    rename("SPECIES" = "INSECT_ID")
  
  ## --------------------------------------------------
  ## Bind and prep annual data
  
  # rbind data from different years
  df <- rbind(df_2021_prepped, df_2022_prepped) %>%
    # change SITE to factor for left_join with HABITAT_CATEGORY created below
    mutate(SITE = as.factor(SITE)) %>%
    # need to convert dates to julian day of year
    mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
    # need to centre and scale the dates as a predictor variable (z-score)
    mutate(SCALED_DATE = center_scale(DATE))
  
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
  
  interval_vector <- as.vector(levels(as.factor(df$INTERVAL)))
  
  visit_vector <- as.vector(levels(as.factor(df$SAMPLING_ROUND_ADJUSTED)))
  
  # find study dimensions from the pooled data
  n_species <- length(species_vector)
  
  n_sites <- length(site_vector)
  
  n_intervals <- length(levels(as.factor(df$INTERVAL)))
  
  n_visits <- length(levels(as.factor(df$SAMPLING_ROUND_ADJUSTED)))
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V
  
  V <- array(data = NA, dim = c(n_species, n_sites, n_intervals, n_visits))
  
  for(i in 1:n_intervals){
    for(j in 1:n_visits){
      
      # iterate this across visits within intervals
      temp <- df_joined %>%
        # filter to indices for interval and visit
        filter(INTERVAL == (i - 1), # intervals start at 0 so need to subtract 1 from i
               SAMPLING_ROUND_ADJUSTED == j) %>% 
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(species_list, by="SPECIES") %>%
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        full_join(site_list, by="SITE") %>%
        separate_rows(SITE, sep = ",") %>%
        # group by SPECIES
        group_by(SPECIES) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(SITE, row, fill = 0) %>%
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(7:24, ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 7:24 represent the columns of each site in the matrix
        #just need the matrix of 0's and 1's
        select(levels(as.factor(df$SITE))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(SPECIES  %in% levels(as.factor(df$SPECIES)))
      

      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      
      # replace NAs for the interval i and visit j with the matrix
      V[1:n_species, 1:n_sites,i,j] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V) <- "numeric"
  
  ## --------------------------------------------------
  ## make an array of the scaled visit dates which will be used as a detection covariate
  
  # need to get all dates from all sites not just those with detections
  # so we need to start from the beginning to pull this information back out
  
  # identify unique visits in 2021
  visits_2021 <- df_2021 %>%
    group_by(SAMPLING_ROUND, SITE) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    mutate(SAMPLING_ROUND_ADJUSTED = SAMPLING_ROUND - 1) %>% # (not using first visits data bc only surveyed for bumblebees)
    mutate(INTERVAL = as.numeric(0)) %>% # year 0
    select(SITE, INTERVAL, SAMPLING_ROUND_ADJUSTED, DATE) 
  
  # identify unique visits in 2022
  visits_2022 <- df_2022 %>%
    group_by(SAMPLING_ROUND, SITE) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    mutate(SAMPLING_ROUND_ADJUSTED = SAMPLING_ROUND) %>% 
    mutate(INTERVAL = as.numeric(1)) %>% # year 1
    filter(SAMPLING_ROUND_ADJUSTED != "7") %>% # currently the model can't handle missing data
    # need to find a way to update the model, otherwise we can't use round 7 where only some sites were visited, 
    # and no sites visited in 2021.
    select(SITE, INTERVAL, SAMPLING_ROUND_ADJUSTED, DATE)
  
  # rbind data from different years
  # and create a scaled date
  visits_all_years <- rbind(visits_2021, visits_2022) %>%
    # change SITE to factor for left_join with HABITAT_CATEGORY created below
    mutate(SITE = as.factor(SITE)) %>%
    filter(!is.na(SAMPLING_ROUND_ADJUSTED)) %>% # there was an NA row for a potential missing specimen
    # need to convert dates to julian day of year
    mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
    # need to centre and scale the dates as a predictor variable (z-score)
    mutate(SCALED_DATE = center_scale(DATE)) 
  
  # initialize a visit array of dimensions n_sites*n_intervals*n_visits
  scaled_date_array <- array(NA, dim =c(n_sites, n_intervals, n_visits))
  
  for(i in 1:n_visits){
    
        temp <- visits_all_years %>%
          filter(SAMPLING_ROUND_ADJUSTED == i) %>%
          select(SITE, INTERVAL, SCALED_DATE) %>%
          pivot_wider(names_from = INTERVAL, values_from = SCALED_DATE)
        
        scaled_date_array[,,i] <- as.matrix(temp[,2:3])
        
  }
  
  
  ## --------------------------------------------------
  ## Return stuff
  return(list(
    
    n_species = n_species, ## number of species
    n_sites = n_sites, ## number of sites
    n_intervals = n_intervals, ## number of occupancy intervals (years)
    n_visits = n_visits, ## number of samples per year
    V = V, # detection array
    species = species_vector, # alphabetized species names
    sites = site_vector, # alphabetized site names
    intervals = interval_vector, # ordered vector of intervals
    visits = visit_vector, # ordered vector of visits
    date_scaled = scaled_date_array, # detection covariate (array of scaled julian date of visit)
    habitat_category = HABITAT_CATEGORY # occupancy and detection covariate (sites mowed or meadow)
    
  ))
  
}
