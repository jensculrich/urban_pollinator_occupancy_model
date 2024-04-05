library(rstan)
library(viridis)
library(gridExtra)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 1 # >=
filter_nonnative_woody = FALSE
my_data <- process_raw_data(min_unique_detections, filter_nonnative_woody)

stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out_binary_habitat.rds")

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

habitat_binary = TRUE # if true, build the herb plots with x axis range from 0, 1

## --------------------------------------------------
### Pull out data

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_years <- my_data$n_years # number of surveys 
n_visits <- my_data$n_visits
n_years_minus1 <- n_years - 1
species <- seq(1, n_species, by=1)
sites <- seq(1, n_sites, by=1)
years <- seq(1, n_years_minus1, by=1)
date_scaled <- my_data$date_scaled
habitat_type <- my_data$habitat_category
species_interaction_metrics <- my_data$species_interaction_metrics
d_scaled <- species_interaction_metrics$d_scaled
d <- species_interaction_metrics$d
species_names <- my_data$species
site_names <- my_data$sites
herbaceous_flowers_scaled <- my_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_data$woody_flowers_scaled
original_herb_abundance <- my_data$original_herb_abundance
original_woody_abundance <- my_data$original_woody_abundance

original_woody_diversity <- my_data$original_woody_diversity

## --------------------------------------------------
### Pull out model output

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## Prediction for processes across range of herbaceous and perennial flower abundance across range of specialization (d')
# Use bin sizes and colours that match the histogram. Will need to convert between scaled and unscaled

#library(RColorBrewer)
#n_bins = 11
#my_palette <- palette(brewer.pal(n = n_bins, name = "RdYlBu"))
#seq <- seq(-2, 5, 0.65)

n_bins = 12
my_palette <- palette(viridis(n = n_bins, option = "C"))

seq_original <- seq(min(species_interaction_metrics$d), 
           max(species_interaction_metrics$d), 
           length.out=n_bins)

seq_scaled <- seq(min(species_interaction_metrics$d_scaled), 
                    max(species_interaction_metrics$d_scaled), 
                    length.out=n_bins)

## --------------------------------------------------
### Make a histogram of species specialization. Colour the bins.

hist(d)
par(mar = c(5, 5, 3, 3)) # Set the margin on all sides to 2
hist(d, axes = TRUE, xlim = c(0, 1),
     ylab = "Number of species", xlab = "Specialization (d')", 
     main = "", 
     col = my_palette, breaks = seq_original, 
     freq=TRUE, right = FALSE,
     cex.axis=1.5, cex.lab=2)

## --------------------------------------------------
# global plotting options

tmp <- as.data.frame(stan_out) # take estimates from each HMC step as a df
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?
pred_length <- 100 # divide each covariate axis into some interval steps


xlabel = "log(herbaceous flower abundance / 20m^2)"

## --------------------------------------------------
# get some prediction secquence data

original_herb <- seq(min(original_herb_abundance), 
                     max(original_herb_abundance), 
                     length.out = pred_length) # reasonable prediction range
herb_pred <- (original_herb - mean(original_herb_abundance)) / sd(original_herb_abundance) # unscale the dates

original_woody <- seq(min(original_woody_abundance), 
                     max(original_woody_abundance), 
                     length.out = pred_length) # reasonable prediction range
woody_pred <- (original_woody - mean(original_woody_abundance)) / sd(original_woody_abundance) # unscale the dates

original_woody <- seq(min(original_woody_diversity), 
                      max(original_woody_diversity), 
                      length.out = pred_length) # reasonable prediction range
woody_pred <- (original_woody - mean(original_woody_diversity)) / sd(original_woody_diversity) # unscale the dates



## --------------------------------------------------
# global plotting options (alternative if using binary hab categories)

if(habitat_binary == TRUE){
  original_herb <- as.numeric(0:1) # reasonable prediction range 
  herb_pred <- original_herb
  xlabel <- ""
}

## --------------------------------------------------
# colonization

## --------------------------------------------------
# fill predictions for the community and specialization bins
# based on the parameter estimates and some set of covariate data

predC <- array(NA, dim=c(pred_length, n_samp, 2)) # community means
predSpec <- array(NA, dim=c(pred_length, n_samp, 2, n_bins)) # trends by specialization bin

for(i in 1:n_samp){
  
  # community means don't depend on specialization bin
  predC[,i,1] <- ilogit( # herbaceous flower trend
    # gamma0 +
    tmp[i,8] + 
      # gamma_herbaceous_flowers*x +
      tmp[i,10]*herb_pred
  )
  predC[,i,2] <- ilogit( # woody flower trend
    # gamma0 +
    tmp[i,8] + 
      
      tmp[i,10]*1 + # for an enhanced park
      
      # gamma_herbaceous_flowers*x +
      tmp[i,11]*woody_pred
  )
  
  # predictions by specialization bin
  for(j in 1:n_bins){
    predSpec[,i,1,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,8] + 
        # gamma_herbaceous_flowers*x +
        tmp[i,10] * herb_pred +
        # gamma_specialization*d[i] + 
        tmp[i,12] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,13] * seq_scaled[j] * herb_pred
    )
    
    predSpec[,i,2,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,8] + 
        
        tmp[i,10]*1 + # for a enhanced park
        
        # gamma_herbaceous_flowers*x +
        tmp[i,11] * woody_pred +
        # gamma_specialization*d[i] + 
        tmp[i,12] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,14] * seq_scaled[j] * woody_pred
    )

  }
}

# posterior means by community average 
criC <- apply(predC, c(1,3), function(x) quantile(x, 
                    prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

# posterior means for specialization specific 
criSpec <- apply(predSpec, c(1,3,4), function(x) quantile(x, prob = c(0.25, 0.5, 0.75)))

#-------------------------------------------------------------------------------

# community plot - gamma herb
herb_df <- as.data.frame(cbind(original_herb, criC[3,,1], 
                               criC[1,,1], criC[5,,1],
                               criC[2,,1], criC[4,,1])) %>%
  rename("gamma_herb_community_mean" = "V2",
         "gamma_herb_community_lower95" = "V3",
         "gamma_herb_community_upper95" = "V4",
         "gamma_herb_community_lower50" = "V5",
         "gamma_herb_community_upper50" = "V6")

p <- ggplot(data = herb_df, aes(original_herb, gamma_herb_community_mean)) +
  geom_ribbon(aes(
    ymin=gamma_herb_community_lower50, 
    ymax=gamma_herb_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=gamma_herb_community_lower95, 
    ymax=gamma_herb_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("colonization rate \n(community average)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
                     ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# specialization bin plot - gamma herb

herb_df_spec <- as.data.frame(cbind(original_herb, 
                                    rep(1:n_bins, each=pred_length),
                                    as.vector(criSpec[2,,1,1:12]), 
                                    as.vector(criSpec[1,,1,1:12]),
                                    as.vector(criSpec[3,,1,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
 # filter(specialization_bin == "12")

p2 <- ggplot(data = herb_df_spec, aes(original_herb, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("colonization rate \n(by diet specialization)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2

grid.arrange(p, p2, ncol=2)

#-------------------------------------------------------------------------------

# community plot - gamma woody
woody_df <- as.data.frame(cbind(original_woody, criC[3,,2], 
                               criC[1,,2], criC[5,,2],
                               criC[2,,2], criC[4,,2])) %>%
  rename("gamma_woody_community_mean" = "V2",
         "gamma_woody_community_lower95" = "V3",
         "gamma_woody_community_upper95" = "V4",
         "gamma_woody_community_lower50" = "V5",
         "gamma_woody_community_upper50" = "V6")

q <- ggplot(data = woody_df, aes(original_woody, gamma_woody_community_mean)) +
  geom_ribbon(aes(
    ymin=gamma_woody_community_lower50, 
    ymax=gamma_woody_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=gamma_woody_community_lower95, 
    ymax=gamma_woody_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("colonization rate \n(community average)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
q

# specialization bin plot - gamma herb

woody_df_spec <- as.data.frame(cbind(original_woody, 
                                    rep(1:n_bins, each=pred_length),
                                    as.vector(criSpec[2,,2,1:12]), 
                                    as.vector(criSpec[1,,2,1:12]),
                                    as.vector(criSpec[3,,2,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
# filter(specialization_bin == "12")

q2 <- ggplot(data = woody_df_spec, aes(original_woody, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("colonization rate \n(by diet specialization)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
q2

grid.arrange(p, p2, q, q2, ncol=2)

## --------------------------------------------------
# persistence

## --------------------------------------------------
# fill predictions for the community and specialization bins
# based on the parameter estimates and some set of covariate data

predC <- array(NA, dim=c(pred_length, n_samp, 2)) # community means
predSpec <- array(NA, dim=c(pred_length, n_samp, 2, n_bins)) # trends by specialization bin

for(i in 1:n_samp){
  
  # community means don't depend on specialization bin
  predC[,i,1] <- ilogit( # herbaceous flower trend
    # gamma0 +
    tmp[i,15] + 
      # gamma_herbaceous_flowers*x +
      tmp[i,17]*herb_pred
  )
  predC[,i,2] <- ilogit( # woody flower trend
    # gamma0 +
    tmp[i,15] + 
      
      tmp[i,17]*1 + # for an enhanced park
      
      # gamma_herbaceous_flowers*x +
      tmp[i,18]*woody_pred
  )
  
  # predictions by specialization bin
  for(j in 1:n_bins){
    predSpec[,i,1,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,15] + 
        # gamma_herbaceous_flowers*x +
        tmp[i,17] * herb_pred +
        # gamma_specialization*d[i] + 
        tmp[i,19] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,20] * seq_scaled[j] * herb_pred
    )
    
    predSpec[,i,2,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,15] + 
        
        tmp[i,17]*1 + # for an enhanced park
        
        # gamma_herbaceous_flowers*x +
        tmp[i,18] * woody_pred +
        # gamma_specialization*d[i] + 
        tmp[i,19] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,21] * seq_scaled[j] * woody_pred
    )
    
  }
}

# posterior means by community average 
criC <- apply(predC, c(1,3), function(x) quantile(x, 
                                                  prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

# posterior means for specialization specific 
criSpec <- apply(predSpec, c(1,3,4), function(x) quantile(x, prob = c(0.25, 0.5, 0.75)))

#-------------------------------------------------------------------------------

# community plot - phi herb
herb_df <- as.data.frame(cbind(original_herb, criC[3,,1], 
                               criC[1,,1], criC[5,,1],
                               criC[2,,1], criC[4,,1])) %>%
  rename("phi_herb_community_mean" = "V2",
         "phi_herb_community_lower95" = "V3",
         "phi_herb_community_upper95" = "V4",
         "phi_herb_community_lower50" = "V5",
         "phi_herb_community_upper50" = "V6")

r <- ggplot(data = herb_df, aes(original_herb, phi_herb_community_mean)) +
  geom_ribbon(aes(
    ymin=phi_herb_community_lower50, 
    ymax=phi_herb_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=phi_herb_community_lower95, 
    ymax=phi_herb_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("persistence rate \n(community average)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
r

# specialization bin plot - phi herb

herb_df_spec <- as.data.frame(cbind(original_herb, 
                                    rep(1:n_bins, each=pred_length),
                                    as.vector(criSpec[2,,1,1:12]), 
                                    as.vector(criSpec[1,,1,1:12]),
                                    as.vector(criSpec[3,,1,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
# filter(specialization_bin == "12")

r2 <- ggplot(data = herb_df_spec, aes(original_herb, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("persistence rate \n(by diet specialization)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
r2

grid.arrange(r, r2, ncol=2)

#-------------------------------------------------------------------------------

# community plot - gamma woody
woody_df <- as.data.frame(cbind(original_woody, criC[3,,2], 
                                criC[1,,2], criC[5,,2],
                                criC[2,,2], criC[4,,2])) %>%
  rename("phi_woody_community_mean" = "V2",
         "phi_woody_community_lower95" = "V3",
         "phi_woody_community_upper95" = "V4",
         "phi_woody_community_lower50" = "V5",
         "phi_woody_community_upper50" = "V6")

s <- ggplot(data = woody_df, aes(original_woody, phi_woody_community_mean)) +
  geom_ribbon(aes(
    ymin=phi_woody_community_lower50, 
    ymax=phi_woody_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=phi_woody_community_lower95, 
    ymax=phi_woody_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("persistence rate \n(community average)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
s

# specialization bin plot - gamma herb

woody_df_spec <- as.data.frame(cbind(original_woody, 
                                     rep(1:n_bins, each=pred_length),
                                     as.vector(criSpec[2,,2,1:12]), 
                                     as.vector(criSpec[1,,2,1:12]),
                                     as.vector(criSpec[3,,2,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
# filter(specialization_bin == "12")

s2 <- ggplot(data = woody_df_spec, aes(original_woody, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("persistence rate \n(by diet specialization)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
s2

grid.arrange(r, r2, s, s2, ncol=2)

## --------------------------------------------------
# initial occupancy

## --------------------------------------------------
# fill predictions for the community and specialization bins
# based on the parameter estimates and some set of covariate data

predC <- array(NA, dim=c(pred_length, n_samp, 2)) # community means
predSpec <- array(NA, dim=c(pred_length, n_samp, 2, n_bins)) # trends by specialization bin

for(i in 1:n_samp){
  
  # community means don't depend on specialization bin
  predC[,i,1] <- ilogit( # herbaceous flower trend
    # gamma0 +
    tmp[i,1] + 
      # gamma_herbaceous_flowers*x +
      tmp[i,3]*herb_pred
  )
  predC[,i,2] <- ilogit( # woody flower trend
    # gamma0 +
    tmp[i,1] + 
      
      tmp[i,3]*1 + # for a meadow park
      
      # gamma_herbaceous_flowers*x +
      tmp[i,4]*woody_pred
  )
  
  # predictions by specialization bin
  for(j in 1:n_bins){
    predSpec[,i,1,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,1] + 
        # gamma_herbaceous_flowers*x +
        tmp[i,3] * herb_pred +
        # gamma_specialization*d[i] + 
        tmp[i,5] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,6] * seq_scaled[j] * herb_pred
    )
    
    predSpec[,i,2,j] <- ilogit( # herbaceous flower trend
      # gamma0 +
      tmp[i,1] + 
        
        tmp[i,3]*1 + # for a meadow park
        
        # gamma_herbaceous_flowers*x +
        tmp[i,4] * woody_pred +
        # gamma_specialization*d[i] + 
        tmp[i,5] * seq_scaled[j] +  
        # gamma_interaction1 *d[i] * x + 
        tmp[i,7] * seq_scaled[j] * woody_pred
    )
    
  }
}

# posterior means by community average 
criC <- apply(predC, c(1,3), function(x) quantile(x, 
                                                  prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

# posterior means for specialization specific 
criSpec <- apply(predSpec, c(1,3,4), function(x) quantile(x, prob = c(0.25, 0.5, 0.75)))

#-------------------------------------------------------------------------------

# community plot - psi1 herb
herb_df <- as.data.frame(cbind(original_herb, criC[3,,1], 
                               criC[1,,1], criC[5,,1],
                               criC[2,,1], criC[4,,1])) %>%
  rename("psi1_herb_community_mean" = "V2",
         "psi1_herb_community_lower95" = "V3",
         "psi1_herb_community_upper95" = "V4",
         "psi1_herb_community_lower50" = "V5",
         "psi1_herb_community_upper50" = "V6")

t <- ggplot(data = herb_df, aes(original_herb, psi1_herb_community_mean)) +
  geom_ribbon(aes(
    ymin=psi1_herb_community_lower50, 
    ymax=psi1_herb_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi1_herb_community_lower95, 
    ymax=psi1_herb_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("initial occurrence rate \n(community average)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
t

# specialization bin plot - phi herb

herb_df_spec <- as.data.frame(cbind(original_herb, 
                                    rep(1:n_bins, each=pred_length),
                                    as.vector(criSpec[2,,1,1:12]), 
                                    as.vector(criSpec[1,,1,1:12]),
                                    as.vector(criSpec[3,,1,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
# filter(specialization_bin == "12")

t2 <- ggplot(data = herb_df_spec, aes(original_herb, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("initial occurrence rate \n(by diet specialization)") +
  xlab(xlabel) +
  scale_x_continuous(limits = c(-0.1,1.1), breaks = c(0,1), labels = c(
    "0" = "conventional", "1" = "enhanced")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
t2

#-------------------------------------------------------------------------------

# community plot - gamma woody
woody_df <- as.data.frame(cbind(original_woody, criC[3,,2], 
                                criC[1,,2], criC[5,,2],
                                criC[2,,2], criC[4,,2])) %>%
  rename("psi1_woody_community_mean" = "V2",
         "psi1_woody_community_lower95" = "V3",
         "psi1_woody_community_upper95" = "V4",
         "psi1_woody_community_lower50" = "V5",
         "psi1_woody_community_upper50" = "V6")

u <- ggplot(data = woody_df, aes(original_woody, psi1_woody_community_mean)) +
  geom_ribbon(aes(
    ymin=psi1_woody_community_lower50, 
    ymax=psi1_woody_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi1_woody_community_lower95, 
    ymax=psi1_woody_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("initial occurrence rate \n(community average)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
u

# specialization bin plot - gamma herb

woody_df_spec <- as.data.frame(cbind(original_woody, 
                                     rep(1:n_bins, each=pred_length),
                                     as.vector(criSpec[2,,2,1:12]), 
                                     as.vector(criSpec[1,,2,1:12]),
                                     as.vector(criSpec[3,,2,1:12]))) %>%
  mutate(specialization_bin = as.factor(V2)) %>%
  rename("mean" = "V3",
         "lower50" = "V4",
         "upper50" = "V5")

#herb_df_spec <- herb_df_spec %>%
# filter(specialization_bin == "12")

u2 <- ggplot(data = woody_df_spec, aes(original_woody, mean, fill=specialization_bin)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.1) +
  geom_line(aes(colour=specialization_bin), size=2, lty=1) +
  xlim(c(min(original_woody), max(original_woody))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("initial occurrence rate \n(by diet specialization)") +
  xlab("log(woody flower abundance / 100m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent 
  ) +
  scale_fill_manual(values=my_palette) +
  scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
u2

grid.arrange(t, t2, u, u2, ncol=2)
