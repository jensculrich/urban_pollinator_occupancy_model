library(reshape2)
library(tidyverse)

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
  filter(unique_detections > 1) %>%
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

# visualize the real detection matrix (number of detections per species/site/year)
V2 <- array(data = 0, dim = c(n_species, n_sites, n_years))

for(i in 1:n_species){
  for(j in 1:n_sites){
    for(k in 1:n_years){
      V2[i,j,k] = sum(V[i,j,k,1:n_visits])
    }
  }
}

df_melted_year1 <- melt(V2[1:n_species,1:n_sites,1])
colnames(df_melted_year1)[1] <- "species"
colnames(df_melted_year1)[2] <- "site"
colnames(df_melted_year1)[3] <- "detections"

df_melted_year2 <- melt(V2[1:n_species,1:n_sites,2])
colnames(df_melted_year2)[1] <- "species"
colnames(df_melted_year2)[2] <- "site"
colnames(df_melted_year2)[3] <- "detections"

df_melted_year3 <- melt(V2[1:n_species,1:n_sites,3])
colnames(df_melted_year3)[1] <- "species"
colnames(df_melted_year3)[2] <- "site"
colnames(df_melted_year3)[3] <- "detections"

detections_heatmap_plot_year1 <-
  ggplot(df_melted_year1, aes(as.factor(site), as.factor(species))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "site number", y = "species number") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Total detections per species/site 2021")

detections_heatmap_plot_year2 <-
  ggplot(df_melted_year2, aes(as.factor(site), as.factor(species))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "site number", y = "species number") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Total detections per species/site 2022")

detections_heatmap_plot_year3 <-
  ggplot(df_melted_year3, aes(as.factor(site), as.factor(species))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "site number", y = "species number") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Total detections per species/site 2021")

gridExtra::grid.arrange(
  detections_heatmap_plot_year1, detections_heatmap_plot_year2, detections_heatmap_plot_year3,
  nrow=1)
