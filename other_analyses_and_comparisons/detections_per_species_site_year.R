library(reshape2)
library(tidyverse)

max_unique_detections = 500 # can lower this to filter off the top if too many species to plot in one page
min_unique_detections = 1

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
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
  filter(!is.na(NO_ID_RESOLUTION))



## --------------------------------------------------
## Transform into detection non detection data

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
  filter(unique_detections < max_unique_detections) %>%
  ungroup() %>%
  
  # change SITE to factor for left_join with HABITAT_CATEGORY created below
  mutate(SITE = as.factor(SITE))

## Get unique species and sites
# create an alphabetized list of all species encountered across all sites*intervals*visits
species_list <- df %>%
  group_by(SPECIES) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SPECIES) # extract species names column as vector

# create an alphabetized list of all sites
site_list <- df %>%
  group_by(SITE) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SITE)

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


# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(df$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites)))

# join HABITAT_CATEGORY of sites with df by SITE
df_joined <- left_join(df, habitats_df, by = "SITE")



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

habitats_df <- habitats_df %>%
  rename("site" = "V3") %>%
  mutate(site = as.integer(site))

clades <- df %>%
  select(CLADE, SPECIES) %>%
  group_by(CLADE, SPECIES) %>%
  slice(1)

species_list <- species_list %>%
  mutate("species" = seq(1:nrow(.))) %>%
  left_join(., clades)

df_melted_year1 <- melt(V2[1:n_species,1:n_sites,1])
colnames(df_melted_year1)[1] <- "species"
colnames(df_melted_year1)[2] <- "site"
colnames(df_melted_year1)[3] <- "detections"
df_melted_year1 <- left_join(df_melted_year1, habitats_df) %>%
  left_join(., species_list) %>%
  mutate(SITE = fct_reorder(SITE, (as.integer(HABITAT_CATEGORY))))

df_melted_year2 <- melt(V2[1:n_species,1:n_sites,2])
colnames(df_melted_year2)[1] <- "species"
colnames(df_melted_year2)[2] <- "site"
colnames(df_melted_year2)[3] <- "detections"
df_melted_year2 <- left_join(df_melted_year2, habitats_df) %>%
  left_join(., species_list) %>%
  mutate(SITE = fct_reorder(SITE, (as.integer(HABITAT_CATEGORY))))

df_melted_year3 <- melt(V2[1:n_species,1:n_sites,3])
colnames(df_melted_year3)[1] <- "species"
colnames(df_melted_year3)[2] <- "site"
colnames(df_melted_year3)[3] <- "detections"
df_melted_year3 <- left_join(df_melted_year3, habitats_df) %>%
  left_join(., species_list) %>%
  mutate(SITE = fct_reorder(SITE, (as.integer(HABITAT_CATEGORY))))

# hline depends on number of species, for min 1 detection right now is 51.5

detections_heatmap_plot_year1 <-
  ggplot(df_melted_year1, aes(SITE, fct_rev((fct_reorder(SPECIES, desc(CLADE)))))) +
                              #fct_reorder(SPECIES, desc(SPECIES)))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "Site", y = "Species") +
  theme_bw() +
  geom_vline(xintercept = 9.5, linewidth = 1, colour="goldenrod", linetype = "longdash") +
  geom_hline(yintercept = 51.5, linewidth = 1, colour="orange3", linetype = "longdash") +
  theme(#legend.position = "bottom",
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size=5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Detections per species/site 2021")

plot(detections_heatmap_plot_year1)

detections_heatmap_plot_year2 <-
  ggplot(df_melted_year2, aes(SITE, fct_rev((fct_reorder(SPECIES, desc(CLADE)))))) +
  #fct_reorder(SPECIES, desc(SPECIES)))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "Site", y = "Species") +
  theme_bw() +
  geom_vline(xintercept = 9.5, linewidth = 1, colour="goldenrod", linetype = "longdash") +
  geom_hline(yintercept = 51.5, linewidth = 1, colour="orange3", linetype = "longdash") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size=5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Detections per species/site 2022")

detections_heatmap_plot_year3 <-
  ggplot(df_melted_year3, aes(SITE, fct_rev((fct_reorder(SPECIES, desc(CLADE)))))) +
  #fct_reorder(SPECIES, desc(SPECIES)))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "Site", y = "Species") +
  theme_bw() +
  geom_vline(xintercept = 9.5, linewidth = 1, colour="goldenrod", linetype = "longdash") +
  geom_hline(yintercept = 51.5, linewidth = 1, colour="orange3", linetype = "longdash") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size=5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Detections per species/site 2023")

gridExtra::grid.arrange(
  detections_heatmap_plot_year1, detections_heatmap_plot_year2, detections_heatmap_plot_year3,
  nrow=1)

