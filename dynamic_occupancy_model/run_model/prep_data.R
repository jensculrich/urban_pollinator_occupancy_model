library(tidyverse)
library(lubridate)

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
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = ifelse(YEAR==1, SAMPLING_ROUND - 1, SAMPLING_ROUND))

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
  filter(unique_detections > 1)

# Species I had to fix by hand: 
# Lasioglossum cressonii; Lasioglossum brunnieventre, Bombus melanopygus, Bombus vosnesenskii, Eumerus sp.
# Merodon equestris, Syrphus opinator

ggplot(mydata_filtered, aes(x=fct_infreq(SPECIES))) +
  geom_bar(stat = "count") +
  labs(x = "Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Total detections per species")
  
ggplot(mydata_filtered_2_plus_unique, aes(x=fct_infreq(SPECIES))) +
  geom_bar(stat = "count") +
  labs(x = "Species") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Total detections per species (for species detected >1 unique visits)")


## --------------------------------------------------
## Transform into detection non detection data

# rbind data from different years
df <- mydata_filtered %>%
  group_by(SPECIES, SITE, YEAR, SAMPLING_ROUND) %>%
  slice(1) %>%
  ungroup() %>%
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

year_vector <- as.vector(levels(as.factor(df$YEAR)))

visit_vector <- as.vector(levels(as.factor(df$SAMPLING_ROUND)))

# find study dimensions from the pooled data
n_species <- length(species_vector)

n_sites <- length(site_vector)

n_years <- length(levels(as.factor(df$YEAR)))

n_visits <- length(levels(as.factor(df$SAMPLING_ROUND)))

## --------------------------------------------------
## Now we are ready to create the detection matrix, V

V <- array(data = NA, dim = c(n_species, n_sites, n_years, n_visits))

for(i in 1:n_years){
  for(j in 1:n_visits){
    
    # iterate this across visits within intervals
    temp <- df_joined %>%
      # filter to indices for interval and visit
      filter(YEAR == i, SAMPLING_ROUND == j) %>% 
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
      mutate_at(15:32, ~replace(., . > 1, 1))%>%
      # if more columns are added these indices above^ might need to change
      # 15:32 represent the columns of each site in the matrix
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
    
    # replace NAs for the year i and visit j with the matrix
    V[1:n_species, 1:n_sites,i,j] <- temp_matrix[1:n_species, 1:n_sites]
    
  }
}

class(V) <- "numeric"
sum(V)
