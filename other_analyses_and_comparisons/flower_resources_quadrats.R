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
mydata <- read.csv("./data/flower_resources_quadrats.csv")

str(mydata)

# perform some initial filters on the unfinished prelim data
mydata <- mydata %>% 
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
  # NA's indicate NO FLOWERS WERE IN BLOOM 
  mutate(NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0))

## Get unique species and sites
# create an alphabetized list of all species encountered across all sites*intervals*visits
species_list <- mydata %>%
  group_by(SPECIES) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SPECIES) # extract species names column as vector

# Species spelling/names I needed to fix:
# white vicia -> Vicia hirsuta, Veronica serpyllfolia -> Veronica serpyllifolia, 
# Trifolium incarnata, Trifolium incarnatatum, Trifolium incanum -> Trifolium incarnatum
# Symphyotrichum douglasii, Senecio vulgaris, Rorippa palustris, Ranunculus acris, 
# Plantago majus, Plantago sp., Nepeta cataria, Mycelis muralis, Melilotus albus,
# Matricaria discoidea, Galium aparine, Erodium cicutarium, Echinacea purpurea, Convolvulus arvensis
# Brassica rapa, Brassica oleracea -> Brassica sp., 
# Bellis perennis, Anaphalis margaritacea, Achillea millefolium

# create an alphabetized list of all sites
site_list <- mydata %>%
  group_by(SITE) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SITE)

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
  scale_x_discrete(labels = c("Control","Enhanced")) 
  # ggtitle("Flower resources in lawn/meadow space")

p  

## --------------------------------------------------
## Model abundance

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

# with random effects and positive only continuous distribution (Gamma) 
temp <- abundance_df_reduced %>%
  mutate(log_total_flower_abundance2 = log_total_flower_abundance + 0.001) # there are some issues with the zeros being "negative"?

(abundance_fit3 <- rstanarm::stan_glmer(log_total_flower_abundance2 ~ HABITAT_CATEGORY + YEAR + 
                                          SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                        family=Gamma(link = "log"),
                                        data=temp)) 

plot(abundance_fit3)

pp_check(abundance_fit3) # posterior pred check by eye looks not very good

posterior_vs_prior(abundance_fit3, pars = c("(Intercept)", "HABITAT_CATEGORY1"))


## --------------------------------------------------
## Let's compare diversity

diversity_df <- mydata %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(species_richness = n_distinct(SPECIES)) %>%
  mutate(sample_id = cur_group_id()) %>% 
  mutate(SITE = as.factor(SITE)) %>%
  ungroup() %>%
  
  mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%

  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(shannon_index = vegan::diversity(log_NUM_FLORAL_UNITS, "shannon")) %>%
  mutate(simpsons_index = vegan::diversity(log_NUM_FLORAL_UNITS, "simpson")) %>%
  ungroup()

# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(diversity_df$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites))) %>%
  rename("SITE_NUMBER" = "V3")

# join HABITAT_CATEGORY of sites with df by SITE
diversity_df_joined <- left_join(diversity_df, habitats_df, by = "SITE")

# for analysis and plotting just want one row per site/year/visit/
diversity_df_reduced <- diversity_df_joined %>%
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  slice(1) %>%
  select(YEAR, SAMPLING_ROUND, SITE, SITE_NUMBER, DATE, 
         HABITAT_CATEGORY, species_richness, shannon_index, simpsons_index) %>%
  ungroup()

q <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=species_richness, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (species richness)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

q

r <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=shannon_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Shannon index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

r

s <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=simpsons_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Simpson's index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

s

gridExtra::grid.arrange(
  p, q, r, s,
  nrow=1)

## --------------------------------------------------
## Model diversity

diversity_df_reduced <- diversity_df_reduced %>%
  mutate(SAMPLING_ROUND = as.factor(SAMPLING_ROUND),
         YEAR = as.factor(YEAR),
         SITE = as.factor(SITE),
         HABITAT_CATEGORY = as.factor(HABITAT_CATEGORY)) %>%
  # need to convert dates to julian day of year
  mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
  # need to centre and scale the dates as a predictor variable (z-score)
  mutate(SCALED_DATE = center_scale(DATE), 
         SCALED_DATE_SQ = SCALED_DATE^2)

## --------------------------------------------------
## Model species richness 

library(lme4)

summary(glm1 <- glm(species_richness ~ HABITAT_CATEGORY, 
             poisson(link = "log"),
             data=diversity_df_reduced))

summary(glm2 <- glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                 SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
               poisson(link = "log"),
               data=diversity_df_reduced))


(richness_fit <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY, 
                                    poisson(link = "log"),
                                    data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(richness_fit)

pp_check(richness_fit)

posterior_vs_prior(richness_fit)


(richness_fit1 <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY, 
                                     neg_binomial_2(link = "log"),
                                     data=diversity_df_reduced))

# posterior predictive checks for poisson look a bit better
plot(richness_fit1)

pp_check(richness_fit1)

posterior_vs_prior(richness_fit1)

# with random effects
(richness_fit2 <- rstanarm::stan_glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                                         SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                       poisson(link = "log"), 
                                       data=diversity_df_reduced)) 

plot(richness_fit2)

pp_check(richness_fit2)

posterior_vs_prior(richness_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

## --------------------------------------------------
## Model species diversity (shannon index)  - gaussian error

(shannon_fit <- rstanarm::stan_glm(shannon_index ~ HABITAT_CATEGORY, 
                                   data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(shannon_fit)

pp_check(shannon_fit)

posterior_vs_prior(shannon_fit)


# with random effects
(shannon_fit2 <- rstanarm::stan_glmer(shannon_index ~ HABITAT_CATEGORY + YEAR + 
                                        SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                      data=diversity_df_reduced)) 

plot(shannon_fit2)

pp_check(shannon_fit2)

posterior_vs_prior(shannon_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))


## -----------------------------------------------------------------------------
## Replicate comparisons with "visited" plant species only

rm(list=ls()) 
gc()

## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## --------------------------------------------------
## Floral Resources in the quadrats from transects of the lawn / meadow space

# read data
mydata <- read.csv("./data/flower_resources_quadrats.csv")

str(mydata)

# perform some initial filters on the unfinished prelim data
mydata <- mydata %>% 
  
  # Reduce sampling rounds in year 1 by 1 (they start at 2 since we did a weird prelim survey first)
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND)))) %>%
  
  # NA's indicate NO FLOWERS WERE IN BLOOM 
  mutate(NUM_FLORAL_UNITS = replace_na(NUM_FLORAL_UNITS, 0))

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

ggplot(plants_visited, aes(x=PLANT_NETTED_FROM_SCI_NAME, y=log_n)) +
  geom_col() +
  ylim(c(0, 7)) +
  labs(x = "Plant species", y="log(number of detected interactions)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ggtitle("Total interactions per species")

mydata_subset <- mydata %>%
  filter(SPECIES %in% plants_visited$PLANT_NETTED_FROM_SCI_NAME)

# create an alphabetized list of all species encountered across all sites*intervals*visits
species_list_reduced <- plants_visited %>%
  select(PLANT_NETTED_FROM_SCI_NAME) # extract species names column as vector


## --------------------------------------------------
## Abundance

abundance_df <- mydata_subset %>%
  
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
  scale_x_discrete(labels = c("Control","Enhanced")) 
# ggtitle("Flower resources in lawn/meadow space")

p  

## --------------------------------------------------
## Model abundance

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

# effect size for flower abundance is about twice as large when we've filetered to "visited" plants

## --------------------------------------------------
## Let's compare diversity

diversity_df <- mydata_subset %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(species_richness = n_distinct(SPECIES)) %>%
  mutate(sample_id = cur_group_id()) %>% 
  mutate(SITE = as.factor(SITE)) %>%
  ungroup() %>%
  
  mutate(log_NUM_FLORAL_UNITS = log(NUM_FLORAL_UNITS + 1)) %>%
  
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(shannon_index = vegan::diversity(log_NUM_FLORAL_UNITS, "shannon")) %>%
  mutate(simpsons_index = vegan::diversity(log_NUM_FLORAL_UNITS, "simpson")) %>%
  ungroup()

# add habitat_category of sites
# make a new dataframe with site names and cats
SITE <- levels(as.factor(diversity_df$SITE))
HABITAT_CATEGORY <- c(0,0,1,1,0,0,1,1,1, # 9 mowed sites 
                      0,1,1,1,0,0,0,1,0) # and 9 meadow sites

habitats_df <- as.data.frame(cbind(SITE, HABITAT_CATEGORY, seq(1:n_sites))) %>%
  rename("SITE_NUMBER" = "V3")

# join HABITAT_CATEGORY of sites with df by SITE
diversity_df_joined <- left_join(diversity_df, habitats_df, by = "SITE")

# for analysis and plotting just want one row per site/year/visit/
diversity_df_reduced <- diversity_df_joined %>%
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  slice(1) %>%
  select(YEAR, SAMPLING_ROUND, SITE, SITE_NUMBER, DATE, 
         HABITAT_CATEGORY, species_richness, shannon_index, simpsons_index) %>%
  ungroup()

q <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=species_richness, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (species richness)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

q

r <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=shannon_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Shannon index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

r

s <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=simpsons_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Simpson's index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

s

gridExtra::grid.arrange(
  p, q, r, s,
  nrow=1)

## --------------------------------------------------
## Model diversity

diversity_df_reduced <- diversity_df_reduced %>%
  mutate(SAMPLING_ROUND = as.factor(SAMPLING_ROUND),
         YEAR = as.factor(YEAR),
         SITE = as.factor(SITE),
         HABITAT_CATEGORY = as.factor(HABITAT_CATEGORY)) %>%
  # need to convert dates to julian day of year
  mutate(DATE = lubridate::yday(as.Date(lubridate::mdy(DATE)))) %>%
  # need to centre and scale the dates as a predictor variable (z-score)
  mutate(SCALED_DATE = center_scale(DATE), 
         SCALED_DATE_SQ = SCALED_DATE^2)

## --------------------------------------------------
## Model species richness 

summary(glm1 <- glm(species_richness ~ HABITAT_CATEGORY, 
                    poisson(link = "log"),
                    data=diversity_df_reduced))

summary(glm2 <- glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                        SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
                      poisson(link = "log"),
                      data=diversity_df_reduced))


(richness_fit <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY, 
                                    poisson(link = "log"),
                                    data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(richness_fit)

pp_check(richness_fit)

posterior_vs_prior(richness_fit)


(richness_fit1 <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY, 
                                     neg_binomial_2(link = "log"),
                                     data=diversity_df_reduced))

# posterior predictive checks for poisson look a bit better
plot(richness_fit1)

pp_check(richness_fit1)

posterior_vs_prior(richness_fit1)

# with random effects
(richness_fit2 <- rstanarm::stan_glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                                         SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                       poisson(link = "log"), 
                                       data=diversity_df_reduced)) 

plot(richness_fit2)

pp_check(richness_fit2)

posterior_vs_prior(richness_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

## --------------------------------------------------
## Model species diversity (shannon index)  - gaussian error

(shannon_fit <- rstanarm::stan_glm(shannon_index ~ HABITAT_CATEGORY, 
                                   data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(shannon_fit)

pp_check(shannon_fit)

posterior_vs_prior(shannon_fit)


# with random effects
(shannon_fit2 <- rstanarm::stan_glmer(shannon_index ~ HABITAT_CATEGORY + YEAR + 
                                        SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                      data=diversity_df_reduced)) 

plot(shannon_fit2)

pp_check(shannon_fit2)

posterior_vs_prior(shannon_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))


#-------------------------------------------------------------------------------
## How do flower resources change across years in the two habitat types?

## --------------------------------------------------
## Abundance
year_names <- list("1" = "2021",
     "2" = "2022",
     "3" = "2023")

year_labeller <- function(variable,value){
  return(year_names[value])
}

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
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) +
  scale_x_discrete(labels = c("Control","Enhanced")) +
  facet_wrap(~YEAR, labeller=year_labeller)
# ggtitle("Flower resources in lawn/meadow space")

p  

## --------------------------------------------------
## Model abundance

summary(glm1 <- glm(log_total_flower_abundance ~ HABITAT_CATEGORY:YEAR, 
                    data=abundance_df_reduced))

summary(glm2 <- lmer(log_total_flower_abundance ~ HABITAT_CATEGORY:YEAR + 
                       SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
                     data=abundance_df_reduced))


(abundance_fit <- rstanarm::stan_glm(log_total_flower_abundance ~ HABITAT_CATEGORY:YEAR, 
                                     data=abundance_df_reduced))

plot(abundance_fit)

pp_check(abundance_fit)

posterior_vs_prior(abundance_fit)


# with random effects
(abundance_fit2 <- rstanarm::stan_glmer(log_total_flower_abundance ~ HABITAT_CATEGORY:YEAR + 
                                          SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                        data=abundance_df_reduced)) 

plot(abundance_fit2)

pp_check(abundance_fit2)

posterior_vs_prior(abundance_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

## --------------------------------------------------
## plot diversity

q <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=species_richness, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (species richness)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) +
  scale_x_discrete(labels = c("Control","Enhanced"))  +
  facet_wrap(~YEAR, labeller=year_labeller)

q

r <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=shannon_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Shannon index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) +
  scale_x_discrete(labels = c("Control","Enhanced"))  +
  facet_wrap(~YEAR, labeller=year_labeller)

r

s <- ggplot(diversity_df_reduced, aes(x=HABITAT_CATEGORY, y=simpsons_index, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (Simpson's index)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) +
  scale_x_discrete(labels = c("Control","Enhanced"))  +
  facet_wrap(~YEAR, labeller=year_labeller)

s

gridExtra::grid.arrange(
  p, q, r, s,
  nrow=1)

## --------------------------------------------------
## model diversity

## --------------------------------------------------
## Model species richness 

summary(glm1 <- glm(species_richness ~ HABITAT_CATEGORY:YEAR, 
                    poisson(link = "log"),
                    data=diversity_df_reduced))

summary(glm2 <- glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                        SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
                      poisson(link = "log"),
                      data=diversity_df_reduced))


(richness_fit <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY:YEAR, 
                                    poisson(link = "log"),
                                    data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(richness_fit)

pp_check(richness_fit)

posterior_vs_prior(richness_fit)


(richness_fit1 <- rstanarm::stan_glm(species_richness ~ HABITAT_CATEGORY:YEAR, 
                                     neg_binomial_2(link = "log"),
                                     data=diversity_df_reduced))

# posterior predictive checks for poisson look a bit better
plot(richness_fit1)

pp_check(richness_fit1)

posterior_vs_prior(richness_fit1)

# with random effects
(richness_fit2 <- rstanarm::stan_glmer(species_richness ~ HABITAT_CATEGORY:YEAR + 
                                         SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                       poisson(link = "log"), 
                                       data=diversity_df_reduced)) 

plot(richness_fit2)

pp_check(richness_fit2)

posterior_vs_prior(richness_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

## --------------------------------------------------
## Model species diversity (shannon index)  - gaussian error

(shannon_fit <- rstanarm::stan_glm(shannon_index ~ HABITAT_CATEGORY:YEAR, 
                                   data=diversity_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(shannon_fit)

pp_check(shannon_fit)

posterior_vs_prior(shannon_fit)


# with random effects
(shannon_fit2 <- rstanarm::stan_glmer(shannon_index ~ HABITAT_CATEGORY:YEAR + 
                                        SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                      data=diversity_df_reduced)) 

plot(shannon_fit2)

pp_check(shannon_fit2)

posterior_vs_prior(shannon_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

#-------------------------------------------------------------------------------
## plot site-specific flower resources

abundance_sites <- abundance_df_reduced %>%
  mutate(SITE = fct_reorder(SITE, as.integer(HABITAT_CATEGORY)))
  
p <- ggplot(abundance_sites, aes(x=SITE, y=log_total_flower_abundance, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "log(flower resource abundance)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) #+
  #scale_x_discrete(labels = c("Control","Enhanced")) +
  #facet_wrap(~YEAR, labeller=year_labeller)

p  

diversity_sites <- diversity_df_reduced %>%
  mutate(SITE = fct_reorder(SITE, as.integer(HABITAT_CATEGORY)))

q <- ggplot(diversity_sites, aes(x=SITE, y=species_richness, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Flower resource diversity (species richness)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust=1)) #+
  #scale_x_discrete(labels = c("Control","Enhanced"))  +
  #facet_wrap(~YEAR, labeller=year_labeller)

q

gridExtra::grid.arrange(
  p, q,
  nrow=1)

#-------------------------------------------------------------------------------
## how do plant identities compare between park types?

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
  arrange(desc(log_total_abundance_all_sites)) %>%
  mutate(row = row_number()) %>%
  select(SPECIES, row)

test <- left_join(test, temp, by = "SPECIES") %>%
  mutate(SPECIES = fct_reorder(SPECIES, (row)))

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
        axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
  facet_wrap(test$HABITAT_CATEGORY, ncol=1, labeller=habitat_labeller)
