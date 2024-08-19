library(tidyverse)
library(rstanarm)
library(ggplot2)

### Pollen limitation analysis
## Started 9.22.22 by J Ulrich

# This file has a self contained analysis, using calling STAN models through rstanarm
# Separately (for practice), I will conduct a parallel analysis using a manually
# written STAN model.  

## SITE_TREATMENT
# Plants were placed in arrays at parks with ("treatment") vs without ("control")
# flower enhancements when in flower. 
# n=5 unique control parks and n=6 unique treatment parks (SITE)

## FLOWER_TREATMENT
# On each plant - to account for variation in resource limitation - one flower was
# designated as ("treatment") and received a pollen supplementation, while another flower
# designated as ("control") and experienced and ambient pollination environment.
# generally we used 4 plants per pot, unless there were not 4 or more plants in the
# pot that produced 2 or more healthy flowers during the one week period where plant arrays
# were installed at sites. 

## SEEDS_PRODUCED
# I let the fruits develop for ~3 weeks after the experiment and then clipped and
# bagged them. Fruits were dissected and number of seeds (the response variable) were 
# counted in the lab.

# I then filter the data set to successfully treated pairs, where both plants 
# were open and receptive during the exposure at the field site, the pollen supplementation
# was applied with confidence that the stigma was receptive, and the supplemented
# flower was neither aborted, split, or destroyed.

# I calculate a pollen limitation index for each untreated (control) flower using 
# using it's paired, treated flower from the same plant with the formula:
# PL_INDEX = 1 - (((treatment - control) / treatment))

# I then attempt to examine the effect of the site treatment on PL_INDEX
# PL_INDEX turns out to be bimodally distributed. Flowers from both site
# types follow this general pattern. I interpret this as meaning that either the
# flower was visited by a viable pollinator and made ~as many seeds as the supplemented flower OR
# was not visited by a viable pollinator and made a few seeds from self-pollination.

# Because the response of PL_INDEX is bimodally distributed, I transform the response
# into a binary outcome: was pollen limited or was not pollen limited.
# I use a threshold of 0.5 (from visual interpretation of the bimodal response)
# with plants < 0.5 being pollen limited (1) or > 0.5 not pollen limited (0).
# I then use a bayesian mixed effects glm (binary family, logit link) to estimate
# the effect of site type (mowed or not mowed) on the probability of the 
# outcome that a flower is pollen limited (1).

## expected outcomes
# A negative effect of no mowing on the pollen limitation outcome would support 
# the ecosystem services spillover hypothesis, that the flower enhancements support
# more insects that end up increase the delivery of the ES of pollination in 
# nearby areas.

###-----------------------------------------------------------------------------
## Analysis with simulated data

# fill with data simulation

###-----------------------------------------------------------------------------
## Analysis with real data

df <- read.csv("./pollen_limitation_experiment/data/clarkia_pollination_data_2022.csv")

###-----------------------------------------------------------------------------
## Basic data summary

# How many pairs did we tag and place in the field
# (Including those we would have tagged but didn't have a flower for in a pot) 
pairs <- 0.5 * ( 
  # multiply by .5, because nrow counts both flowers from a plant
  df %>%
    nrow())

# How many pairs flowered in the window and were successfully treated
# with supplemental pollen during the window
successfully_treated_pairs <- 0.5 * ( 
  # multiply by .5, because nrow counts both flowers from a plant
  df %>%
  filter(BOTH_FLOWERS_OPEN_AND_POLLEN_APPLIED == "Y") %>%
  nrow())

# remove flowers with uncounted seed pods
successfully_treated_pairs_recovered <- 0.5 * ( 
  # multiply by .5, because nrow counts both flowers from a plant
  df %>%
    filter(BOTH_FLOWERS_OPEN_AND_POLLEN_APPLIED == "Y",
           BOTH_CAPSULES_RECOVERED == "Y") %>%
    nrow())

# number pairs where either a) didn't both flower during the experiment or 
# b) both flowered but the stigma of treatment flower was not receptive when we visited the site
(pairs - successfully_treated_pairs)


# number pairs where we didn't recover both capsules either due to
# a) capsule split during the ripening phase (maybe from the heat?)
# b) capsule was removed, presumably by herbivory or weather (capsule fully missing)
# c) labelling on the seed packet was unclear (two packets with the same label)
# d) seed packet was not found in the seed packet bin (assistants misplaced?)
(successfully_treated_pairs - successfully_treated_pairs_recovered)

###-----------------------------------------------------------------------------
## Data prep

max_PL_accepted <- 2

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
  filter(PL_INDEX < max_PL_accepted) %>%
  
  # Now we also likely won't be able to accept true 0 values if we log link the data
  # For now, these will be subbed with infitisemely small - but postive - PL_INDEX values
  mutate(PL_INDEX = ifelse(PL_INDEX == 0, 0.001, PL_INDEX))

###-----------------------------------------------------------------------------
## Visualize prepped data (And transform into a binary response)

# We can see that we have a bimodal response
# Consequently we will likely have to recode the variable into..
# PL_INDEX < threshold (pollen supplementation made a difference) 
# or
# PL_INDEX > threshold (pollen supplementation has no effect)
PL_response_threshold <- 0.5

# Plot density of response
(ggplot(df_prepped) +
   
   geom_density(aes(x=PL_INDEX, fill=SITE_TREATMENT), alpha=0.5) +
   scale_x_continuous(name = "1 - Pollen Limitation",
                      limits = c(0, 2), 
                      breaks = c(0, .5, 1, 1.5, 2),
                      label = c("0", "0.5", "1", "1.5", "2")) +
   scale_y_continuous(name = "Density", limits = c(0,1),
                      breaks = c(0, .5, 1),
                      label = c("0", "0.5", "1")) +
   scale_fill_discrete(name = "", 
                       breaks = c("control", "treatment"),
                       labels = c("control", "restored")) + 
   theme_bw() +
   theme(legend.text=element_text(size=14),
         axis.text.x = element_text(size = 14),
         axis.text.y =element_text(size = 14),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
   geom_vline(xintercept=PL_response_threshold, lty="dashed")
 
)

# What is the mean PL_INDEX (given the filters we've already imposed on failed flowers)
mean(df_prepped$PL_INDEX)



# Here 1 indicates pollen limited
# and 0 indicates no pollen limitation 
df_binary <- df_prepped %>%
  mutate(PL_INDEX_BINARY = ifelse(PL_INDEX < PL_response_threshold, 1, 0))

# Plot density of response
(ggplot(df_binary) +
    
    geom_density(aes(x=PL_INDEX_BINARY, fill=SITE_TREATMENT), alpha=0.5)
  
)

df_binary_plot_data <- df_binary %>%
  group_by(SITE_TREATMENT) %>%
  add_tally(PL_INDEX_BINARY == 1) %>%
  add_tally(PL_INDEX_BINARY == 0) %>%
  slice(1) %>%
  mutate(PROPORTION_POLLEN_LIMITED = n / (n+nn)) %>%
  mutate(PROPORTION_NOT_POLLEN_LIMITED = nn / (n+nn)) %>%
  select(SITE_TREATMENT, 
         PROPORTION_POLLEN_LIMITED, PROPORTION_NOT_POLLEN_LIMITED) %>%
  pivot_longer(!SITE_TREATMENT, names_to = "type", values_to = "proportion")

# Plot density of response
(ggplot(df_binary_plot_data, aes(x=type, y=proportion, fill=SITE_TREATMENT)) +
    
    geom_bar(stat="identity", width=.5, position = "dodge") +
    scale_x_discrete(name="", breaks = c("PROPORTION_NOT_POLLEN_LIMITED", "PROPORTION_POLLEN_LIMITED"),
                       labels=c("Not pollen limited", "Pollen limited")) +
    scale_y_continuous(name="Proportion of flowers", 
                       limits = (c(0,1)), breaks = c(0, .25, .5, .75, 1),
                       label = c("0%", "25%", "50%", "75%", "100%")) +
    scale_fill_discrete(name = "Site type", 
                        breaks = c("control", "treatment"),
                        labels = c("mowed", "unmowed")) +
    theme_bw() +
    theme(legend.text=element_text(size=10),
        axis.text.x = element_text(size = 12),
        axis.text.y =element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
  
)  

###-----------------------------------------------------------------------------
## Fit model

# run glmer with gamma distribution
stan_fit <- stan_glmer(PL_INDEX_BINARY ~ 
                              SITE_TREATMENT + 
                              (1|SITE),
                            family = binomial(link = "logit"),
                            data = df_binary) 

# confirm sufficient mixing of the chains
(trace <- plot(stan_fit, "trace", pars = c("(Intercept)",
                                           "SITE_TREATMENTtreatment")))
# confirm stable exploration of the possible parameter space
(pairs(stan_fit, pars = c("(Intercept)",
                                         "SITE_TREATMENTtreatment")))

saveRDS(stan_fit, "./stan_fit.RDS")
 
###-----------------------------------------------------------------------------
## View model outputs

# values less than zero decrease the odds that a plant is pollen limited
print(stan_fit)  
round(posterior_interval(stan_fit, prob = 0.95), 2)

# posterior distribution for the parameters describing the uncertainty 
# related to unknown parameter values:
plot(stan_fit)
pplot <- plot(stan_fit, pars = c("(Intercept)",
                                 "SITE_TREATMENTtreatment",
                                 "Sigma[SITE:(Intercept),(Intercept)]")
              )
pplot + 
  geom_vline(xintercept = 0) +
  scale_y_discrete(breaks=c("(Intercept)",
                            "SITE_TREATMENTtreatment",
                            "Sigma[SITE:(Intercept),(Intercept)]"
                            ),
                   labels=c("Intercept",
                            "Effect of no mow",
                            "Sigma (site)"
                   )) +
  scale_x_continuous(
    name="Posterior Parameter Estimate \n (Values > 0 increase odds of pollen limitation; \n Values < 0 decrease odds of pollen limitation)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y =element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

plot_title <- ggplot2::ggtitle("Posterior Distributions")
plot(stan_fit, "dens_overlay", pars = c("(Intercept)",
                                        "SITE_TREATMENTtreatment")) + 
  plot_title

###-----------------------------------------------------------------------------
## LOO comparison

# leave one out cross validation
(loo1 <- loo(stan_fit, save_psis = TRUE))

stan_fit0 <- update(stan_fit, formula = PL_INDEX_BINARY ~ (1|SITE), QR = FALSE, refresh=0)
(loo0 <- loo(stan_fit0))
loo_compare(loo0, loo1)
# the baseline model is (marginally) better.
# covariate (site treatment) DOES NOT contain clearly useful information for predictions.

###------------------------------------------------------------------------------
## Posterior Predictive Check

# can the model can reproduce the mean observed in the real data?
# (where the mean is the proportion of plants that are 
# pollen limited accourding to our pollen limitation threshold)
(test <- pp_check(stan_fit, plotfun = "stat", binwidth = 0.01))

# The histogram shows the distributions of the fraction of simulated data 
# that resulted in pollen limitation in each of the posterior draws.

# can the model can reproduce the patterns observed in the real data
# where we have most plants being pollinated relatively sufficiently
# and about 1/3 being pollen limited?
pp_check(stan_fit, nreps=100)