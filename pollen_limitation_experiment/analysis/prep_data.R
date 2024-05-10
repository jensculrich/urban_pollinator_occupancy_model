library(tidyverse)

## Prepare real data for mixed-effects, logistic regression
#  started Dec. 14, 2022, J Ulrich

# Format binary outcome data where.. 
# the outcome is: 
# 0 - a flower is not pollen limited, OR
# 1 - a flower is pollen limited
# A random intercept for planter pot, nested in site will account for heterogeneity
# in nested sample groups.

prep_data <- function(max_PL_accepted,
                      PL_response_threshold){
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  ###-----------------------------------------------------------------------------
  ## Read real data
  
  df <- read.csv("./pollen_limitation_experiment/data/clarkia_pollination_data_2022.csv")
  
  ###-----------------------------------------------------------------------------
  ## Data prep
  
  # prepare data set for analysis
  df_prepped <- df %>%
    
    # remove flowers with from pairs not treated or not recovered
    filter(BOTH_FLOWERS_OPEN_AND_POLLEN_APPLIED == "Y",
           BOTH_CAPSULES_RECOVERED == "Y") %>%
    
    # I added excel formulas to view PL and PL INDEX, but let's drop and recalculate here
    # to make sure that the formula is clear and reproducible.
    # We will only use the PL_INDEX, which I use to represent the number of seeds
    # produced by a flower under ambient pollination conditions, relative to 
    # to the number of seeds produced by a flower on the same plant that receives a complete pollen load
    # PL_INDEX = 1 - ((Treatment - Control) / Control)
    # Should always be > 0
    # In theory, PL_INDEXshould be less than 1, but may be greater than 1 due to 
    # differences in ovules per capsule or treatment failure.
    select(-PL, -PL_INDEX) %>% # drop the values calculated in excel
    select(SITE_TREATMENT, SITE, POT_NUMBER, FLOWER_TREATMENT, SEEDS_PRODUCED) %>% # drop unneeded columns
    # add a plant ID variable
    mutate(PLANT_ID = rep(1:(.5*nrow(.)), each = 2)) %>%
    # now spread by values for seeds produced within Plant ID
    pivot_wider(names_from = FLOWER_TREATMENT, values_from = SEEDS_PRODUCED) %>%
    # and calculate PL_INDEX
    mutate(PL_INDEX = 1 - ((treatment - control) / treatment)) %>%
    
    # We will remove plants where the treatment produced 0 seeds,
    # presumably this is due to either a failure of the treatment or 
    # or some growth deformity/stress of the plant that completely prevented capsule from developing
    # should be clear here about how removing these plants will effect the results
    # but as is, cannot produce a PL_INDEX from these plants because cannot divide by 0 to get the ratio
    filter(treatment != 0) %>%
    
    # We will also remove plants from a maximum acceptabe PL_INDEX
    # Plants with an index above this value made significantly more seeds for untreated (ambient)
    # flowers than on supplemented flowers. This could be due to overapplication of pollen
    # failure to apply the pollen when the stigma was actually receptive,
    # or some growth deformity/stress of the plant that prevented capsule from developing
    # Again, remember to be transparent about this value if we choose to retain it,
    # and how altering the threshold might affect the results.
    filter(PL_INDEX < max_PL_accepted)
  
  ###-----------------------------------------------------------------------------
  ## Transform from a continuous to binary outcome
  
  # Here 1 indicates pollen limited
  # and 0 indicates no pollen limitation 
  df_binary <- df_prepped %>%
    mutate(PL_INDEX_BINARY = ifelse(PL_INDEX < PL_response_threshold, 1, 0))
  
  ###-----------------------------------------------------------------------------
  ## Extract additional info from the df
  
  # siteLookup (which level 3 cluster is level 2 cluster grouped in?)
  siteLookup <- df_binary %>%
    group_by(as.factor(POT_NUMBER)) %>%
    slice(1) %>%
    pull(SITE)
  
  SITES <- df_binary %>%
    group_by(SITE) %>%
    slice(1) %>%
    pull(SITE)
  
  # site treatment covariate 
  x <- df_binary %>%
    # treatment as 1, control as 0
    mutate(SITE_TREATMENT = (as.numeric(as.factor(SITE_TREATMENT)) - 1)) %>%
    pull(SITE_TREATMENT) 
  
  # site treatment covariate 
  y <- df_binary %>%
    pull(PL_INDEX_BINARY) 
  
  ## --------------------------------------------------
  ## Get pollinator data (so that we can then filter flower abundance down to the "pollinator plants")
  
  # read data
  mydata <- read.csv("./pollen_limitation_experiment/data/pollinator_data.csv")
  
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
  
  ## --------------------------------------------------
  ## Get site X year flowering plant abundance data (herbaceous plants in quadrats)
  
  # read data
  plant_data <- read.csv("./pollen_limitation_experiment/data/flower_resources_quadrats.csv")
  
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
  
  ## now calculate plant abundance as for the regular data set BUT
  # we only want from year 2022 and the surveys before and after the pollination experiment
  # and only from the subset of sites that we were able to do the pollination experiment at
  
  # abundance
  plant_abundance_df <- plant_data_subset %>%
    
    filter(YEAR == 2) %>% # get year of the pollination experiment
    filter(SAMPLING_ROUND %in% c(4, 5)) %>% # filter to mid summer (we surveyed before and after the experiment)
    filter(SITE %in% SITES) %>% # we only did the experiment at 11/18 sites
    
    # calculate abundance per site visit
    group_by(SITE, SAMPLING_ROUND) %>%
    mutate(total_flower_abundance = sum(NUM_FLORAL_UNITS),
           log_total_flower_abundance = log(total_flower_abundance + 1)) %>%
    group_by(SITE) %>% 
    mutate(avg_flower_abundance = mean(total_flower_abundance),
           avg_log_flower_abundance = mean(log_total_flower_abundance)) %>%
    
    # take one row per site visit
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(herb_abundance_scaled = center_scale(avg_log_flower_abundance)) %>%
    
    # select the necessary stuff
    select(SITE, herb_abundance_scaled, avg_log_flower_abundance)
  
  original_herb_abundance <- plant_abundance_df$avg_log_flower_abundance
  
  df_binary <- left_join(df_binary, plant_abundance_df)
  
  ## --------------------------------------------------
  ## Get site X year flowering plant abundance data (woody plants within survey area)
  
  # read data
  woody_plant_data <- read.csv("./pollen_limitation_experiment/data/flower_resources_woody.csv")
  
  # perform some initial filters on the unfinished prelim data
  woody_plant_data <- woody_plant_data %>% 
    
    # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
    mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
    
    # NA's indicate NO FLOWERS WERE IN BLOOM 
    mutate(SPECIES = replace_na(SPECIES, "No woody flowers"),
           NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0)) %>%
    
    filter(SHRUB_OR_TREE == "y") %>%
    filter(SITE %in% SITES) %>%
    filter(YEAR == 2) %>%
    filter(SAMPLING_ROUND %in% c(4,5)) 
  
  n_surveys = 2
  year = 2
  site_visits <- as.data.frame(cbind(
    rep(SITES, each=1*n_surveys), 
    rep(c(4,5), times=length(SITES)*1), 
    rep(year, each = n_surveys, times=1), 
    rep("no flowering species", n_surveys*length(SITES)), 
    rep(0, n_surveys*length(SITES)))) %>%
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
    filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME) %>%
    select(YEAR, SAMPLING_ROUND, SITE, SPECIES, NUM_FLORAL_UNITS)
  
  # join all possible site visits back in to add zeros
  woody_plant_data_subset <- full_join(woody_plant_data_subset, site_visits,
                                       by = c("YEAR", "SAMPLING_ROUND", "SITE")) 
  
  
  # get mean per site X year
  woody_plant_abundance_df <- woody_plant_data_subset %>%
    mutate(NUM_FLORAL_UNITS.x = replace_na(NUM_FLORAL_UNITS.x, 0),
      log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS.x + 1)) %>%
    group_by(SITE, SAMPLING_ROUND) %>%
    mutate(log_NUM_FLORAL_UNITS_sum = sum(log_NUM_FLORAL_UNITS)) %>%
    group_by(SITE) %>%
    mutate(mean_woody_plant_abundance = mean(log_NUM_FLORAL_UNITS_sum)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # scale the variable
    mutate(woody_plant_abundance_scaled = center_scale(mean_woody_plant_abundance)) %>%
    
    # select the necessary stuff
    select(SITE, woody_plant_abundance_scaled, mean_woody_plant_abundance)
  
  original_woody_abundance <- woody_plant_abundance_df$mean_woody_plant_abundance
  
  df_binary <- left_join(df_binary, woody_plant_abundance_df)
  
  ###-----------------------------------------------------------------------------
  ## Return stuff
  
  return(list(
    
    N = nrow(df_binary), # number of pairs
    n_pots = length(unique(df_binary$POT_NUMBER)),
    pots = df_binary$POT_NUMBER, # vector of pot names
    n_sites = length(unique(df_binary$SITE)), # number of sites
    
    # should add a numeric site name and return both
    sites = df_binary$SITE, # vector of site names
    sites_numeric = as.numeric(as.factor(df_binary$SITE)), # numeric site names
    # need a site name for each pot (length = n_pots)
    siteLookup = siteLookup,
    
    x = x, # site type covariate data
    y = y, # outcome data (binary PL index),
    
    herb_abundance_scaled = df_binary$herb_abundance_scaled,
    woody_abundance_scaled = df_binary$woody_plant_abundance_scaled,
    original_herb_abundance = original_herb_abundance,
    original_woody_abundance = original_woody_abundance,
    
    avg_log_flower_abundance = df_binary$avg_log_flower_abundance
    
  ))
  
}