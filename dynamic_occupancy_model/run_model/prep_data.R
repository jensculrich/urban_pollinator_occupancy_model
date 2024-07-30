library(tidyverse) # organization tools
library(lubridate) # prep survey date data 
library(bipartite) # calculate species interaction metrics
library(rstatix) # conduct a t-test for woody plants

process_raw_data <- function(min_unique_detections, filter_nonnative_woody) {

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
  
  # total detections (species detected _ or more unique occasions)
  ggplot(temp, aes(x=fct_infreq(SPECIES))) +
    geom_bar(stat = "count") +
    labs(x = "", y = "Total detections") +
    theme(axis.text.x = element_text(size = 11, angle = 65, vjust = 1, hjust=1.1),
          axis.text.y = element_text(size = 11)) +
    ggtitle(paste0("Total detections per species (for species detected ", min_unique_detections, " or more unique visits)"))
  
  # unique detections (species detected _ or more unique occasions)
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
  
  clade_list <- df_joined %>%
    group_by(SPECIES) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(CLADE) # extract species names column as vector
  
  # create an alphabetized list of all sites
  site_list <- df_joined %>%
    group_by(SITE) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(SITE) # extract species names column as vector
  
  # get vectors of species, sites, intervals, and visits 
  species_vector <- species_list %>%
    pull(SPECIES)
  
  clade_vector <- clade_list %>%
    pull(CLADE)
  
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
    filter(SAMPLING_ROUND < 7) %>%
    group_by(SAMPLING_ROUND, YEAR, SITE) %>% # one unique row per site*species*visit combination
    slice(1) %>%
    ungroup() %>%
    select(SITE, YEAR, SAMPLING_ROUND, DATE) %>%
    
    # need to convert dates to julian day of year
    mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE))),
           # the first year had an unusually warm spring and things got started earlier
           # but we didn't visit earlier (in fact later) than other years
           # potentially account for this shift in phenology by upping the calendar day of year
           # for year 1?
           DATE_ADJUSTED = ifelse(YEAR==1, DATE+10, DATE)) %>%
    # need to centre and scale the dates as a predictor variable (z-score)
    mutate(SCALED_DATE = center_scale(DATE),
           SCALED_DATE_ADJUSTED = center_scale(DATE_ADJUSTED))
  
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
  
  # scaled date adjusted array
  # initialize a visit array of dimensions n_sites*n_years*n_visits
  scaled_date_adjusted_array <- array(NA, dim =c(n_sites, n_years, n_visits))
  
  for(i in 1:n_visits){
    
    temp <- visits %>%
      filter(SAMPLING_ROUND == i) %>%
      select(SITE, YEAR, SCALED_DATE_ADJUSTED) %>%
      pivot_wider(names_from = YEAR, values_from = SCALED_DATE_ADJUSTED)
    
    if(i %in% 1:6){ # if there were six visits in the year:
      scaled_date_adjusted_array[,,i] <- as.matrix(temp[,2:(n_years + 1)])
      
    } else{ # we didnt survey all sites 7 times in all years so we need to handle differently
      
      temp2 <- left_join(site_list, temp) 
      
      # stan can't handle NA's
      # instead we will pick an unfeasible number so that we don't get confused
      # and then later tell stan to just pass over these values (do not evaluate likelihood)
      temp2[is.na(temp2)] <- -99
      
      temp3 <- pull(temp2[,2])
      
      scaled_date_adjusted_array[,2,i] <- temp3
      scaled_date_adjusted_array[,1,i] <- rep(-99, times=n_sites)
      scaled_date_adjusted_array[,3,i] <- rep(-99, times=n_sites)
    }
    
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
  ## Get site X year flowering plant abundance data (herbaceous plants in quadrats)
  
  # read data
  plant_data <- read.csv("./data/flower_resources_herb_quadrats.csv")
  
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
  
  # abundance
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
  
  plant_abundance_df_wide <- plant_abundance_df %>%
    select(SITE, YEAR, plant_abundance_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = plant_abundance_scaled)
  
  herbaceous_flowers_scaled <- as.matrix(plant_abundance_df_wide[2:(n_years+1)])
  
  # diversity
  plant_diversity_df <- plant_data_subset %>%
    
    # calculate abundance per site visit
    group_by(SITE, YEAR, SAMPLING_ROUND) %>%
    add_tally() %>%
    
    # take one row per site visit
    slice(1) %>%
    ungroup() %>%
    
    # there could be multiple ways to group and calculate plant data
    # here grouping by average abundance per year
    group_by(SITE, YEAR) %>%
    mutate(mean_annual_plant_diversity = mean(n)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(plant_diversity_scaled = center_scale(mean_annual_plant_diversity))
  
  original_herb_diversity <- plant_diversity_df$mean_annual_plant_diversity
    
  plant_diversity_df_wide <- plant_diversity_df %>%
    select(SITE, YEAR, plant_diversity_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = plant_diversity_scaled)
  
  herbaceous_diversity_scaled <- as.matrix(plant_diversity_df_wide[2:(n_years+1)])
  
  # corr between abundance and diversity
  cor.test(plant_diversity_df$mean_annual_plant_diversity, plant_abundance_df$mean_annual_plant_abundance)
  
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
  woody_plant_abundance_df <- woody_plant_data_subset %>%
    mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%
    group_by(SITE, YEAR) %>%
    mutate(mean_annual_woody_plant_abundance = mean(log_NUM_FLORAL_UNITS)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(woody_plant_abundance_scaled = center_scale(mean_annual_woody_plant_abundance)) 
  
  original_woody_abundance <- woody_plant_abundance_df$mean_annual_woody_plant_abundance
  
  woody_plant_data_subset_wide <- woody_plant_abundance_df %>%
    select(SITE, YEAR, woody_plant_abundance_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = woody_plant_abundance_scaled)
  
  woody_flowers_scaled <- as.matrix(woody_plant_data_subset_wide[2:(n_years+1)])
  
  woody_plant_data_subset_wide_original_scale <- woody_plant_abundance_df %>%
    select(SITE, YEAR, mean_annual_woody_plant_abundance) %>%
    pivot_wider(names_from = YEAR, values_from = mean_annual_woody_plant_abundance)
  
  woody_flowers_original_scale <- as.matrix(woody_plant_data_subset_wide_original_scale[2:(n_years+1)])
  
  woody_flowers_scaled_all_years <- rowMeans(woody_flowers_original_scale)
  
  # diversity
  woody_plant_diversity_df_hits <- woody_plant_data_subset %>%
    
    filter(SPECIES != "no flowering species") %>%
    # calculate abundance per site visit
    group_by(SITE, YEAR) %>%
    add_tally() %>%
 
    # take one row per site visit
    slice(1) %>%
    ungroup() %>%
    
    rename("annual_woody_plant_diversity" = "n") 
  
  possible_woody_surveys <- as.data.frame(
    cbind(rep(year_vector, times=n_sites), 
    rep(SITE, each = n_years))) %>%
    rename("YEAR" = "V1",
            "SITE" = "V2") %>%
    mutate(YEAR = as.integer(YEAR))
  
  woody_plant_diversity_df <- full_join(
    woody_plant_diversity_df_hits, possible_woody_surveys) %>%
    mutate(annual_woody_plant_diversity = 
             replace_na(annual_woody_plant_diversity, 0)) %>%

    # scale the log-transformed variable
    mutate(plant_diversity_scaled = center_scale(log(annual_woody_plant_diversity+1))) 
  
  # now sort in similar fashion as before
  woody_plant_diversity_df <- with(woody_plant_diversity_df, 
                                   woody_plant_diversity_df[order(
                                     SITE, YEAR) , ])
  
  original_woody_diversity <- woody_plant_diversity_df$annual_woody_plant_diversity
  
  woody_plant_diversity_df_wide <- woody_plant_diversity_df %>%
    select(SITE, YEAR, plant_diversity_scaled) %>%
    pivot_wider(names_from = YEAR, values_from = plant_diversity_scaled)
  
  woody_diversity_scaled <- as.matrix(woody_plant_diversity_df_wide[2:(n_years+1)])
  
  # corr between abundance and diversity
  cor.test(woody_plant_diversity_df$annual_woody_plant_diversity, woody_plant_abundance_df$mean_annual_woody_plant_abundance)
  
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
  ## t-test for woody plants for control versus herb sites
  
  cat_woody <- as.data.frame(
    cbind(HABITAT_CATEGORY, woody_flowers_scaled_all_years)) %>%
    mutate(HABITAT_CATEGORY = gsub('0', 'control', HABITAT_CATEGORY),
           HABITAT_CATEGORY = gsub('1', 'herb. enhancement', HABITAT_CATEGORY))
  
  (summary_stats <- cat_woody %>%
    group_by(HABITAT_CATEGORY) %>%
    get_summary_stats(woody_flowers_scaled_all_years, type = "mean_sd"))
  
  stat.test <- cat_woody %>% 
    t_test(woody_flowers_scaled_all_years ~ HABITAT_CATEGORY) %>%
    add_significance()
  stat.test
  
  p <- ggplot(cat_woody, aes(as.factor(HABITAT_CATEGORY), woody_flowers_scaled_all_years) )
  p + geom_boxplot() + geom_point() +
    xlab("") +
    ylab("log(woody flower abundance)\naveraged across all visits") +
    annotate("text", x = 1, y = -0.5, size=6, label = paste0("mean = ", summary_stats[1,4])) +
    annotate("text", x = 2, y = -0.5, size=6, label = paste0("mean = ", summary_stats[2,4])) +
    annotate("text", x = 1.5, y = 5, size=6, label = paste0("t = ", signif(pull(stat.test[6]), 5))) +
    annotate("text", x = 1.5, y = 4.5, size=6, label = paste0("p-value = ", stat.test[8])) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size = 18),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
  
  
  ## --------------------------------------------------
  ## Return stuff
  return(list(
    
    n_species = n_species, ## number of species
    n_sites = n_sites, ## number of sites
    n_years = n_years, ## number of occupancy intervals (years)
    n_visits = n_visits, ## number of samples per year
    V = V, # detection array
    species = species_vector, # alphabetized species names
    clade = clade_vector, # anthophila versus syrphidae in same order as species names
    sites = site_vector, # alphabetized site names
    years = year_vector, # ordered vector of intervals
    visits = visit_vector, # ordered vector of visits
    date_scaled = scaled_date_array, # detection covariate (array of scaled julian date of visit)
    date_adjusted_scaled = scaled_date_adjusted_array, 
    habitat_category = HABITAT_CATEGORY, # occupancy and detection covariate (sites mowed or meadow)
    species_interaction_metrics = species_interaction_metrics,
    herbaceous_flowers_scaled = herbaceous_flowers_scaled,
    woody_flowers_scaled = woody_flowers_scaled,
    herbaceous_diversity_scaled = herbaceous_diversity_scaled,
    woody_diversity_scaled = woody_diversity_scaled,
    herbaceous_flowers_by_survey = herbaceous_flowers_by_survey,
    woody_flowers_by_survey = woody_flowers_by_survey,
    flowers_any_by_survey = flowers_any_by_survey,
    
    original_herb_abundance = original_herb_abundance,
    original_woody_abundance = original_woody_abundance,
    original_herb_diversity = original_herb_diversity,
    original_woody_diversity = original_woody_diversity
    
  ))
  
}
