library(tidyverse)
library(lubridate)
library(lme4)
library(rstanarm)
library(vegan)

# will compare floral resources in:
# quadrats from transects of the lawn / meadow space
# offered by perennial shrubs and trees located within the sampling area

# might revisit after filtering by "plants that pollinators ever were recorded interacting with"
# might also add in "1"'s to indicate presence of plants that we netted species off of, but 
# that never occurred within our sampling quadrats (false negatives)

## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## --------------------------------------------------
## Floral Resources in the quadrats from transects of the lawn / meadow space

# read data
mydata <- read.csv("./data/flower_resources_woody.csv")

str(mydata)

# perform some initial filters on the unfinished prelim data
mydata <- mydata %>% 
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
  # NA's indicate NO FLOWERS WERE IN BLOOM 
  mutate(SPECIES = replace_na(SPECIES, "No woody flowers"),
         NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0)) %>%
         #, 
         #SHRUB_OR_TREE = 
           #ifelse(SPECIES == "No woody flowers", "y", SHRUB_OR_TREE)) %>%

  filter(SHRUB_OR_TREE == "y") 

counts <- mydata %>%
  group_by(SPECIES) %>%
  mutate(total = sum(NUM_FLORAL_UNITS)) %>%
  slice(1)


## Get unique species and sites
# create an alphabetized list of all species encountered across all sites*intervals*visits
species_list <- mydata %>%
  group_by(SPECIES) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SPECIES) # extract species names column as vector

# create an alphabetized list of all sites
site_list <- mydata %>%
  group_by(SITE) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SITE)

site_vector <- levels(as.factor(site_list$SITE))

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

# join all possible site visits back in to add zeros
mydata <- full_join(mydata, site_visits) 

nrow(test2 <- mydata %>% 
  group_by(SITE, SAMPLING_ROUND, YEAR) %>%
  slice(1)) ; nrow(site_visits)

## --------------------------------------------------
## Abundance

abundance_df <- mydata %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(total_flower_abundance = sum(NUM_FLORAL_UNITS),
         log_total_flower_abundance = log(total_flower_abundance + 1)) %>%
  # change SITE to factor for left_join with HABITAT_CATEGORY created below
  mutate(SITE = as.factor(SITE))

# find study dimensions from the pooled data
n_species <- length(levels(as.factor(abundance_df$SPECIES)))

n_sites <- length(levels(as.factor(abundance_df$SITE)))

n_years <- length(levels(as.factor(abundance_df$YEAR)))

n_visits <- length(levels(as.factor(abundance_df$SAMPLING_ROUND)))

# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(abundance_df$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites))) %>%
  rename("SITE_NUMBER" = "V3")

# join HABITAT_CATEGORY of sites with df by SITE
abundance_df_joined <- left_join(abundance_df, habitats_df, by = "SITE")

# for analysis and plotting just want one row per site/year/visit/
abundance_df_reduced <- abundance_df_joined %>%
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  slice(1) %>%
  select(YEAR, SAMPLING_ROUND, SITE, SITE_NUMBER, DATE, HABITAT_CATEGORY, total_flower_abundance, log_total_flower_abundance) %>%
  ungroup()

p <- ggplot(abundance_df_reduced, aes(x=HABITAT_CATEGORY, y=log_total_flower_abundance, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "log(flower resource abundance)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) +
  ggtitle("Woody plant flower resources per site visit (any woody plant)")

p  

#-------------------------------------------------------------------------------
## plot with only "visited" plant species

# read pollinator data
pollinator_data <- read.csv("./data/pollinator_data.csv")

plants_visited <- pollinator_data %>%
  group_by(PLANT_NETTED_FROM_SCI_NAME) %>%
  add_tally() %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(PLANT_NETTED_FROM_SCI_NAME, n) %>% 
  mutate(log_n = log(n)) %>%
  filter(PLANT_NETTED_FROM_SCI_NAME != "") %>%
  mutate(PLANT_NETTED_FROM_SCI_NAME = fct_reorder(PLANT_NETTED_FROM_SCI_NAME, desc(n))) %>%
  filter(n >= 5)

mydata_subset <- mydata %>%
  filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME)

# join all possible site visits back in to add zeros
mydata_subset <- full_join(mydata_subset, site_visits) 

nrow(mydata_subset <- mydata_subset %>% 
       group_by(SITE, SAMPLING_ROUND, YEAR) %>%
       slice(1)) ; nrow(site_visits)

abundance_df <- mydata_subset %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(total_flower_abundance = sum(NUM_FLORAL_UNITS),
         log_total_flower_abundance = log(total_flower_abundance + 1)) %>%
  # change SITE to factor for left_join with HABITAT_CATEGORY created below
  mutate(SITE = as.factor(SITE))

# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(mydata$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites))) %>%
  rename("SITE_NUMBER" = "V3")

# join HABITAT_CATEGORY of sites with df by SITE
abundance_df_joined <- left_join(abundance_df, habitats_df, by = "SITE")

# for analysis and plotting just want one row per site/year/visit/
abundance_df_reduced <- abundance_df_joined %>%
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  slice(1) %>%
  select(YEAR, SAMPLING_ROUND, SITE, SITE_NUMBER, DATE, HABITAT_CATEGORY, total_flower_abundance, log_total_flower_abundance) %>%
  ungroup()

p <- ggplot(abundance_df_reduced, aes(x=HABITAT_CATEGORY, y=log_total_flower_abundance, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "log(flower resource abundance)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) +
  ggtitle("Woody plant flower resources per site visit")

p  

#-------------------------------------------------------------------------------
## plot with average per site per year

abundance_df <- mydata_subset %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR) %>%
  mutate(total_flower_abundance = sum(NUM_FLORAL_UNITS),
         log_total_flower_abundance = log(total_flower_abundance + 1)) %>%
  # change SITE to factor for left_join with HABITAT_CATEGORY created below
  mutate(SITE = as.factor(SITE))

# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(mydata$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites))) %>%
  rename("SITE_NUMBER" = "V3")

# join HABITAT_CATEGORY of sites with df by SITE
abundance_df_joined <- left_join(abundance_df, habitats_df, by = "SITE")

# for analysis and plotting just want one row per site/year/visit/
abundance_df_reduced <- abundance_df_joined %>%
  group_by(SITE, YEAR) %>%
  slice(1) %>%
  select(YEAR, SITE, SITE_NUMBER, DATE, HABITAT_CATEGORY, total_flower_abundance, log_total_flower_abundance) %>%
  ungroup()

q <- ggplot(abundance_df_reduced, aes(x=HABITAT_CATEGORY, y=log_total_flower_abundance, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "log(flower resource abundance)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) +
  ggtitle("Average woody plant flower resources per site/year")

q 

gridExtra::grid.arrange(
  p, q, 
  nrow=1)

#-------------------------------------------------------------------------------
## model differences in abundance

# need to join dates back in to do this.

abundance_df_reduced <- abundance_df_reduced %>%
  mutate(SAMPLING_ROUND = as.factor(SAMPLING_ROUND),
         YEAR = as.factor(YEAR),
         SITE = as.factor(SITE),
         HABITAT_CATEGORY = as.factor(HABITAT_CATEGORY)) %>%
  # need to convert dates to julian day of year
  mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
  # need to centre and scale the dates as a predictor variable (z-score)
  mutate(SCALED_DATE = center_scale(DATE), 
         SCALED_DATE_SQ = SCALED_DATE^2)



summary(glm1 <- glm(log_total_flower_abundance ~ HABITAT_CATEGORY, 
                    data=abundance_df_reduced))

summary(glm2 <- lmer(log_total_flower_abundance ~ HABITAT_CATEGORY + YEAR + 
                       SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
                     data=abundance_df_reduced))


(abundance_fit <- rstanarm::stan_glm(log_total_flower_abundance ~ HABITAT_CATEGORY, 
                                     data=abundance_df_reduced))

plot(abundance_fit)

pp_check(abundance_fit)

posterior_vs_prior(abundance_fit)


# with random effects
(abundance_fit2 <- rstanarm::stan_glmer(log_total_flower_abundance ~ HABITAT_CATEGORY + YEAR + 
                                          SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                        data=abundance_df_reduced)) 

plot(abundance_fit2)

pp_check(abundance_fit2)

posterior_vs_prior(abundance_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

#-------------------------------------------------------------------------------
## show variation across sites

#-------------------------------------------------------------------------------
## compare diversity

#-------------------------------------------------------------------------------
## again show variation across sites

#-------------------------------------------------------------------------------
## what kind of woody shrubs were in bloom?
## how much were each of these visited?

test <- abundance_df_joined %>%
  group_by(HABITAT_CATEGORY, SPECIES) %>%
  mutate(total_abundance_all_sites = sum(NUM_FLORAL_UNITS), 
         log_total_abundance_all_sites = log(total_abundance_all_sites)) %>%
  slice(1) %>%
  ungroup() %>%
  #mutate(SPECIES = fct_reorder(SPECIES, desc(log_total_abundance_all_sites))) %>%
  mutate(HABITAT_CATEGORY = as.factor(HABITAT_CATEGORY)) %>%
  
  # add column to indicate whether species occurs in both site types or not
  group_by(SPECIES) %>%
  add_tally() %>%
  ungroup() %>%
  mutate(n = as.factor(n)) 

temp <- abundance_df_joined %>%
  group_by(HABITAT_CATEGORY, SPECIES) %>%
  mutate(total_abundance_all_sites = sum(NUM_FLORAL_UNITS), 
         log_total_abundance_all_sites = log(total_abundance_all_sites)) %>%
  slice(1) %>%
  ungroup() %>% 
  filter(HABITAT_CATEGORY == 1) %>%
  mutate(row = row_number()) %>%
  select(SPECIES, row) 

test <- left_join(test, temp, by = "SPECIES") %>%
  filter(SPECIES != "no flowering species") 
  

str(test)

habitat_names <- list("0" = "control", "1" = "enhanced")

habitat_labeller <- function(variable,value){
  return(habitat_names[value])
}

ggplot(test, aes(x=SPECIES, y=log_total_abundance_all_sites, fill=n)) +
  geom_col() +
  labs(x = "Plant species", y="log(abundance)") +
  scale_fill_brewer(palette="Greens") + theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 60, hjust=1, size = 12),
        strip.text.x = element_text(size = 14)) +
  facet_wrap(test$HABITAT_CATEGORY, ncol=1, labeller=habitat_labeller)
