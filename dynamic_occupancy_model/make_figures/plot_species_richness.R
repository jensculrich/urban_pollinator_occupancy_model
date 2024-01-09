library(rstan)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 1 # >=
my_data <- process_raw_data(min_unique_detections)

stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest


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
d <- species_interaction_metrics$d_scaled
species_names <- my_data$species
site_names <- my_data$sites
herbaceous_flowers_scaled <- my_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_data$woody_flowers_scaled


## --------------------------------------------------
### Pull out model output

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## Use colour palette that is consistent with palette used for other figures

library(viridis)
n_bins = 15
my_palette <- palette(viridis(n = n_bins, option = "C"))
my_palette_reduced <- rev(c(my_palette[2], my_palette[7]))

## --------------------------------------------------
## categorical habitat species richness plot

## --------------------------------------------------
# set up for plot
# parameter means
X <- seq(1:n_years) # number of years

# mean of eco params
Y <- c(
  fit_summary$summary[93,1], # mean species richness year 1
  fit_summary$summary[94,1], # mean species richness year 2
  fit_summary$summary[95,1] # mean species richness year 3
)

# confidence intervals
lower_95 <- c(
  fit_summary$summary[93,4], # mean species richness year 1
  fit_summary$summary[94,4], # mean species richness year 2
  fit_summary$summary[95,4] # mean species richness year 3
)

upper_95 <- c(
  fit_summary$summary[93,8], # mean species richness year 1
  fit_summary$summary[94,8], # mean species richness year 2
  fit_summary$summary[95,8] # mean species richness year 3
)

lower_50 <- c(
  fit_summary$summary[93,5], # mean species richness year 1
  fit_summary$summary[94,5], # mean species richness year 2
  fit_summary$summary[95,5] # mean species richness year 3
)

upper_50 <- c(
  fit_summary$summary[93,7], # mean species richness year 1
  fit_summary$summary[94,7], # mean species richness year 2
  fit_summary$summary[95,7] # mean species richness year 3
)

df_estimates_enhanced <- as.data.frame(cbind(X, Y, 
                                             lower_95, upper_95,
                                             lower_50, upper_50)) %>%
  mutate(habitat = "enhanced")

# mean of eco params
Y <- c(
  fit_summary$summary[90,1], # mean species richness year 1
  fit_summary$summary[91,1], # mean species richness year 2
  fit_summary$summary[92,1] # mean species richness year 3
)

# confidence intervals
lower_95 <- c(
  fit_summary$summary[90,4], # mean species richness year 1
  fit_summary$summary[91,4], # mean species richness year 2
  fit_summary$summary[92,4] # mean species richness year 3
)

upper_95 <- c(
  fit_summary$summary[90,8], # mean species richness year 1
  fit_summary$summary[91,8], # mean species richness year 2
  fit_summary$summary[92,8] # mean species richness year 3
)

lower_50 <- c(
  fit_summary$summary[90,5], # mean species richness year 1
  fit_summary$summary[91,5], # mean species richness year 2
  fit_summary$summary[92,5] # mean species richness year 3
)

upper_50 <- c(
  fit_summary$summary[90,7], # mean species richness year 1
  fit_summary$summary[91,7], # mean species richness year 2
  fit_summary$summary[92,7] # mean species richness year 3
)

df_estimates_control <- as.data.frame(cbind(X, Y, 
                               lower_95, upper_95,
                               lower_50, upper_50)) %>%
  mutate(habitat = "control")

df_estimates <- rbind(df_estimates_enhanced, df_estimates_control)

df_estimates$X <- as.factor(df_estimates$X)

## --------------------------------------------------
## Draw species richness plot

(p <- ggplot(data = df_estimates, aes(X, Y, color = habitat)) +
  scale_x_discrete(name="", breaks = c(1, 2, 3),
                    labels=c("2021", "2022", "2023"
                    )) +
  geom_point(size = 8, position=position_dodge(width=0.5)) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.1,
    position=position_dodge(width=0.5),
    size=1) +
   geom_errorbar(
     aes(ymin = lower_50, ymax = upper_50),
     width = 0,
     position=position_dodge(width=0.5),
     size=2.5) +
  ylim(c(40, 85)) +
  theme_bw() +
  ylab("species richness") +
  scale_color_manual(values=my_palette_reduced) +
  theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)


