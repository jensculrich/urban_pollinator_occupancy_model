library(tidyverse)
library(lubridate)

# read data
mydata <- read.csv("./data/pollinator_data.csv")

## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

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
  mutate(SAMPLING_ROUND = as.integer(ifelse(YEAR==1, as.integer(SAMPLING_ROUND) - 1, as.integer(SAMPLING_ROUND))))%>%
  
  filter(!is.na(NO_ID_RESOLUTION))

## --------------------------------------------------
## Let's compare abundances

abundance_df <- mydata_filtered %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  add_tally() %>%
  rename("abundance" = "n") %>%
  
  # change SITE to factor for left_join with HABITAT_CATEGORY created below
  mutate(SITE = as.factor(SITE))

## Get unique species and sites
# create an alphabetized list of all species encountered across all sites*intervals*visits
species_list <- abundance_df %>%
  group_by(SPECIES) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SPECIES) # extract species names column as vector

# create an alphabetized list of all sites
site_list <- abundance_df %>%
  group_by(SITE) %>%
  slice(1) %>% # take one row per species (the name of each species)
  ungroup() %>%
  select(SITE)


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
  select(YEAR, SAMPLING_ROUND, SITE, SITE_NUMBER, DATE, HABITAT_CATEGORY, abundance) %>%
  ungroup()

p <- ggplot(abundance_df_reduced, aes(x=HABITAT_CATEGORY, y=abundance, fill=HABITAT_CATEGORY)) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=23, size=6) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  labs(title="",x="", y = "Pollinator abundance") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_x_discrete(labels = c("Control","Enhanced")) 

  #geom_hline(yintercept = 10) +
  #geom_hline(yintercept = 14)

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

library(lme4)

(glm1 <- glm(abundance ~ HABITAT_CATEGORY, 
              poisson(link = "log"),
              data=abundance_df_reduced))

(glm2 <- glmer(abundance ~ HABITAT_CATEGORY + YEAR + 
                 SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
               poisson(link = "log"),
               data=abundance_df_reduced))

library(rstanarm)

(abundance_fit <- rstanarm::stan_glm(abundance ~ HABITAT_CATEGORY, 
                                      poisson(link = "log"),
                                      data=abundance_df_reduced))

# the posterior predictive checks for poisson look REALLY BAD 
plot(abundance_fit)

pp_check(abundance_fit)

posterior_vs_prior(abundance_fit)


(abundance_fit1 <- rstanarm::stan_glm(abundance ~ HABITAT_CATEGORY, 
                                      neg_binomial_2(link = "log"),
                                      data=abundance_df_reduced))

# posterior predictive checks for negative binomial look much better
plot(abundance_fit1)

pp_check(abundance_fit1)

posterior_vs_prior(abundance_fit1)

# with random effects
(abundance_fit3 <- rstanarm::stan_glmer(abundance ~ HABITAT_CATEGORY + YEAR + 
                                          SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                       neg_binomial_2(link = "log"), 
                                       data=abundance_df_reduced)) 

plot(abundance_fit3)

pp_check(abundance_fit3)

posterior_vs_prior(abundance_fit3, pars = c("(Intercept)", "HABITAT_CATEGORY1"))

## --------------------------------------------------
## Let's compare diversity

library(vegan)

diversity_df <- mydata_filtered %>%
  
  # calculate abundance per site/visit
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(species_richness = n_distinct(SPECIES)) %>%
  mutate(sample_id = cur_group_id()) %>% 
  mutate(SITE = as.factor(SITE)) %>%
  ungroup() %>%
  
  group_by(SITE, YEAR, SAMPLING_ROUND, SPECIES) %>%
  add_tally() %>%
  ungroup() %>%
  
  group_by(SITE, YEAR, SAMPLING_ROUND) %>%
  mutate(shannon_index = vegan::diversity(n, "shannon")) %>%
  mutate(simpsons_index = vegan::diversity(n, "simpson")) 

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
  labs(title="",x="", y = "Pollinator diversity (species richness)") +
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
  labs(title="",x="", y = "Pollinator diversity (Shannon index)") +
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
  labs(title="",x="", y = "Pollinator diversity (Simpson's index)") +
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

(glm1 <- glm(species_richness ~ HABITAT_CATEGORY, 
             poisson(link = "log"),
             data=diversity_df_reduced))

(glm2 <- glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                 SCALED_DATE + SCALED_DATE_SQ + (1|SITE), 
               poisson(link = "log"),
               data=diversity_df_reduced))

library(rstanarm)

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

# posterior predictive checks for negative binomial look much better
plot(richness_fit1)

pp_check(richness_fit1)

posterior_vs_prior(richness_fit1)

# with random effects
(richness_fit2 <- rstanarm::stan_glmer(species_richness ~ HABITAT_CATEGORY + YEAR + 
                                          SCALED_DATE + SCALED_DATE_SQ + (1|SITE),
                                        neg_binomial_2(link = "log"), 
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


## --------------------------------------------------
## Model species diversity (simpsons index)- beta regression

temp <- diversity_df_reduced %>%
  mutate(simpsons_index2 = simpsons_index + 0.001) # there are some issues with the zeros being "negative"?

(simpsons_fit <- rstanarm::stan_betareg(simpsons_index2 ~ HABITAT_CATEGORY,
                                   data=temp))

# the posterior predictive checks for poisson look REALLY BAD 
plot(simpsons_fit)

pp_check(simpsons_fit) # NOT VERY GOOD

posterior_vs_prior(simpsons_fit)


# how to add random effects to betareg?
(simpsons_fit2 <- rstanarm::stan_betareg(simpsons_index2 ~ HABITAT_CATEGORY + YEAR + 
                                        SCALED_DATE + SCALED_DATE_SQ,
                                      data=temp)) 

plot(simpsons_fit2)

pp_check(simpsons_fit2)

posterior_vs_prior(simpsons_fit2, pars = c("(Intercept)", "HABITAT_CATEGORY1"))