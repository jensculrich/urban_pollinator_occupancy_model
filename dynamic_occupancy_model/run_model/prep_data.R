library(tidyverse)
library(lubridate)

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
    filter(!SPECIES %in% c("Apis mellifera","undetermined", "undetermined/unconfirmed ID"))  %>%
    
    # Filter by SPECIES
    # Remove some others until resolved for final analyses
    filter(!SPECIES %in% c("Bombus sp.","Eupeodes sp.", "Sphaerophoria sp."))  %>%
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND))))
  
  # number of unique species
  (length(unique(mydata_filtered$SPECIES)))
  
  # easy table summary of total detections per species
  detections_species <- mydata_filtered %>%
    group_by(SPECIES) %>%
    count() %>%
    rename(total_detections = n) 
  
  # remove species detected fewer than 2 unique occassions (site/year/visit) for figure
  mydata_filtered_2_plus_unique <- mydata_filtered %>%  
    group_by(SPECIES, SITE, YEAR, SAMPLING_ROUND) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(SPECIES) %>%
    add_tally() %>%
    rename(unique_detections = n) %>%
    filter(unique_detections >= min_unique_detections)
  
  species_2_or_more <- unique(mydata_filtered_2_plus_unique$SPECIES)
  temp <- mydata_filtered %>%
    filter(SPECIES %in% species_2_or_more)
  
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
  ggplot(mydata_filtered_2_plus_unique, aes(x=fct_infreq(SPECIES))) +
    geom_bar(stat = "count") +
    labs(x = "", y = "Unique detections") +
    theme(axis.text.x = element_text(size = 11, angle = 65, vjust = 1, hjust=1.1),
          axis.text.y = element_text(size = 11)) +
    ggtitle(paste0("Unique site/year/visit detections per species (for species detected ", min_unique_detections, " or more unique visits)"))
  
  
  
  ## --------------------------------------------------
  ## Transform into detection non detection data
  
  # rbind data from different years
  df <- mydata_filtered %>%
    
    # convert to binary detections
    group_by(SPECIES, SITE, YEAR, SAMPLING_ROUND) %>%
    slice(1) %>%
    ungroup() %>%
    
    # remove species detected on < 2 unique occassions
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
    
    scaled_date_array[,,i] <- as.matrix(temp[,2:(n_years + 1)])
    
  }
  
  if(length(which(is.na(scaled_date_array))) != 0){
    print("WARNING: missing survey dates for 1 or more site/year/visit combinations!")
  }
  
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
    date_scaled = scaled_date_array, # detection covariate (array of scaled julian date of visit)
    habitat_category = HABITAT_CATEGORY # occupancy and detection covariate (sites mowed or meadow)
    
  ))
  
}
