library(rstan)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 2 # >=
my_data <- process_raw_data(min_unique_detections)

stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out2.rds")

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

reverse_center_scale <- function(x) {
  mean(x) + z * sd(x)
  #(x - mean(x)) / sd(x)
}

## --------------------------------------------------
## Prediction for processes across range of herbaceous and perennial flower abundance across range of specialization (d')
# Use bin sizes and colours that match the histogram. Will need to convert between scaled and unscaled

#library(RColorBrewer)
#n_bins = 11
#my_palette <- palette(brewer.pal(n = n_bins, name = "RdYlBu"))
#seq <- seq(-2, 5, 0.65)

library(viridis)
n_bins = 16
my_palette <- palette(viridis(n = n_bins, option = "C"))
#seq <- seq(-2, 5, 0.45)
seq <- seq(-2.5, 3.5, 0.4)

## --------------------------------------------------
### Make a histogram of species specialization. Colour the bins.

hist(d)
# scaled values
hist(d, axes = TRUE, xlim = c(-2.5, 3.5), 
     ylab = "Frequency", xlab = "Species specialization (Bluthgen's d) (scaled)", 
     main = "", 
     col = my_palette, breaks = seq, 
     freq=TRUE, right = FALSE)

# make one with actual values of d
# unscale
hist(reverse_center_scale(d), axes = TRUE, xlim = c(-2.5, 3.5), 
     ylab = "Frequency", xlab = "Species specialization (Bluthgen's d) ", 
     main = "", 
     col = my_palette, breaks = seq, freq=TRUE, right = FALSE)

## --------------------------------------------------
# colonization

# model: 
# gamma[i,j,k] =
# gamma0 +
# gamma_species[species[i]] + // a species specific intercept
# specialization
# gamma_herbaceous_flowers[species[i]] * herbaceous_flowers_scaled[j,k] + 
# gamma_woody_flowers[species[i]] * woody_flowers_scaled[j,k] 
# interaction 1
# interaction 2

## --------------------------------------------------
# response to herbaceous flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(herbaceous flower abundance) (scaled)",
     ylab = "Pr(colonization)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[8,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[10,1]*x + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[11,1]*0 + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[12,1]*seq[i] +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[13,1] * seq[i] * x) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[14,1] * seq[i] * 0) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
    # gamma0 +
    fit_summary$summary[8,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[10,1]* x + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[11,1]*0 + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[12,1]*0 +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[13,1] * 0 * x) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[14,1] * 0 * 0) 
  ),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# response to woody flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(woody plant flower abundance) (scaled)",
     ylab = "Pr(colonization)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[8,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[10,1]*0 + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[11,1]*x + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[12,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[13,1] * seq[i] * 0) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[14,1] * seq[i] * x) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[8,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[10,1]*0 + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[11,1]*x + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[12,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[13,1] * 0 * 0) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[14,1] * 0 * x) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# persistence

# model: 
# phi[i,j,k] =
# phi0 +
# phi_species[species[i]] + // a species specific intercept
# specia
# phi_herbaceous_flowers[species[i]] * herbaceous_flowers_scaled[j,k] + 
# phi_woody_flowers[species[i]] * woody_flowers_scaled[j,k] 
# interaction1
# interaction2

## --------------------------------------------------
# response to herbaceous flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(herbaceous flower abundance) (scaled)",
     ylab = "Pr(persistence)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[15,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[17,1]*x + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[18,1]*0 + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[19,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[20,1] * seq[i] * x) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[21,1] * seq[i] * 0) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[15,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[17,1]* x + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[18,1]*0 + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[19,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[20,1] * 0 * x) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[21,1] * 0 * 0) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# response to woody flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(woody plant flower abundance) (scaled)",
     ylab = "Pr(persistence)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[15,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[17,1]*0 + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[18,1]*x + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[19,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[20,1] * seq[i] * 0) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[21,1] * seq[i] * x) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[15,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[17,1]*0 + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[18,1]*x + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[19,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[20,1] * 0 * 0) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[21,1] * 0 * x) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# initial occupancy

## --------------------------------------------------
# response to herbaceous flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(herbaceous flower abundance) (scaled)",
     ylab = "Pr(occupancy in year 1)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[1,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[3,1]*x + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[4,1]*0 + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[5,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[6,1] * seq[i] * x) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[7,1] * seq[i] * 0) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[1,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[3,1]* x + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[4,1]*0 + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[5,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[7,1] * 0 * x) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[8,1] * 0 * 0) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# response to woody flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(woody plant flower abundance) (scaled)",
     ylab = "Pr(occupancy in year 1)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # gamma0 +
    fit_summary$summary[1,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[3,1]*0 + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[4,1]*x + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[5,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[6,1] * seq[i] * 0) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[7,1] * seq[i] * x) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[1,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[3,1]*0 + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[4,1]*x + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[5,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[7,1] * 0 * 0) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[8,1] * 0 * x) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# detection
# need to fix this still

## --------------------------------------------------
# response to herbaceous flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(survey-specific flower abundance) (scaled and averaged)",
     ylab = "Pr(detection)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # p0 +
    fit_summary$summary[1,1] + 
      # gamma_herbaceous_flowers*x +
      fit_summary$summary[3,1]*x + 
      # gamma_woody_flowers*x + # set this to mean of zero
      fit_summary$summary[4,1]*0 + 
      # gamma_specialization*d[i] + 
      fit_summary$summary[5,1]*seq[i] +  
      # gamma_interaction1 *d[i] * x + 
      (fit_summary$summary[6,1] * seq[i] * x) +   
      # gamma_interaction2 *d[i] * x + # set this to mean of zero   
      (fit_summary$summary[7,1] * seq[i] * 0) 
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[1,1] + 
    # gamma_herbaceous_flowers*x +
    fit_summary$summary[3,1]* x + 
    # gamma_woody_flowers*x + # set this to mean of zero
    fit_summary$summary[4,1]*0 + 
    # gamma_specialization*d[i] + 
    fit_summary$summary[5,1]*0 +  
    # gamma_interaction1 *d[i] * x + 
    (fit_summary$summary[7,1] * 0 * x) +   
    # gamma_interaction2 *d[i] * x + # set this to mean of zero   
    (fit_summary$summary[8,1] * 0 * 0) 
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")
