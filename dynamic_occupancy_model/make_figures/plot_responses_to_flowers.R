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
## Prediction for processes across range of herbaceous and perennial flower abundance across range of specialization (d')
# Use bin sizes and colours that match the histogram. Will need to convert between scaled and unscaled

#library(RColorBrewer)
#n_bins = 11
#my_palette <- palette(brewer.pal(n = n_bins, name = "RdYlBu"))
#seq <- seq(-2, 5, 0.65)

library(viridis)
n_bins = 16
my_palette <- palette(viridis(n = n_bins, option = "C"))
seq <- seq(-2, 5, 0.45)

## --------------------------------------------------
### Make a histogram of species specialization. Colour the bins.

hist(d)
# scaled values
hist(d, axes = TRUE, xlim = c(-2, 5), ylab = "Frequency", xlab = "Species specialization (Bluthgen's d) (scaled)", main = "", 
     col = my_palette, breaks = seq, freq=TRUE, right = FALSE)
# make one with actual values of d
# unscale
hist(d, axes = TRUE, xlim = c(-2, 5), ylab = "Frequency", xlab = "Species specialization (Bluthgen's d) (scaled)", main = "", 
     col = my_palette, breaks = seq, freq=TRUE, right = FALSE)

## --------------------------------------------------
# colonization

# model: 
# gamma[i,j,k] =
# gamma0 +
# gamma_species[species[i]] + // a species specific intercept
# gamma_herbaceous_flowers[species[i]] * herbaceous_flowers_scaled[j,k] + 
# gamma_woody_flowers[species[i]] * woody_flowers_scaled[j,k] 
# where community mean in gamma_species, gamma_herbaceous_flowers, and gamma_woody_flowers defined as:
# mu_gamma0[i] = delta1_gamma0*d[i]; // baseline colonization rate (centered on 0)
# mu_gamma_herbaceous_flowers[i] = delta0_gamma_herbaceous + delta1_gamma_herbaceous*d[i]; // effect of specialization on effect of habitat on persistence rate
# mu_gamma_woody_flowers[i] = delta0_gamma_woody + delta1_gamma_woody*d[i]; // effect of specialization on effect of habitat on persistence rate


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
    fit_summary$summary[10,1] + 
    # delta1_gamma0*d[i] +
    fit_summary$summary[11,1]*seq[i] + 
    # (delta0_gamma_herbaceous + delta1_gamma_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[13,1] + fit_summary$summary[14,1]*seq[i]) * x) +
    # hold at constant mean woody # (delta0_gamma_woody + delta1_gamma_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[15,1] + fit_summary$summary[16,1]*seq[i]) * 0)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[10,1] + 
    # delta1_gamma0*d[i] +
    fit_summary$summary[11,1]*0 + 
    # (delta0_gamma_herbaceous + delta1_gamma_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[13,1] + fit_summary$summary[14,1]*0) * x) +
    # hold at constant mean woody # (delta0_gamma_woody + delta1_gamma_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[15,1] + fit_summary$summary[16,1]*0) * 0)
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
    fit_summary$summary[10,1] + 
      # delta1_gamma0*d[i] +
      fit_summary$summary[11,1]*seq[i] + 
      # (delta0_gamma_herbaceous + delta1_gamma_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
      ((fit_summary$summary[13,1] + fit_summary$summary[14,1]*seq[i]) * 0) +
      # hold at constant mean woody # (delta0_gamma_woody + delta1_gamma_woody*d[i]) * woody_flowers_scaled[] <- 0 
      ((fit_summary$summary[15,1] + fit_summary$summary[16,1]*seq[i]) * x)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # gamma0 +
  fit_summary$summary[10,1] + 
    # delta1_gamma0*d[i] +
    fit_summary$summary[11,1]*0 + 
    # (delta0_gamma_herbaceous + delta1_gamma_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[13,1] + fit_summary$summary[14,1]*0) * 0) +
    # hold at constant mean woody # (delta0_gamma_woody + delta1_gamma_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[15,1] + fit_summary$summary[16,1]*0) * x)
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

## --------------------------------------------------
# persistence

# model: 
# phi[i,j,k] =
# phi0 +
# phi_species[species[i]] + // a species specific intercept
# phi_herbaceous_flowers[species[i]] * herbaceous_flowers_scaled[j,k] + 
# phi_woody_flowers[species[i]] * woody_flowers_scaled[j,k] 
# where community mean in phi_species, phi_herbaceous_flowers, and phi_woody_flowers defined as:
# mu_phi0[i] = delta1_gamma0*d[i]; // baseline persistence rate (centered on 0)
# mu_phi_herbaceous_flowers[i] = delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]; // effect of specialization on effect of habitat on persistence rate
# mu_phi_woody_flowers[i] = delta0_phi_woody + delta1_phi_woody*d[i]; // effect of specialization on effect of habitat on persistence rate


## --------------------------------------------------
# response to herbaceous flower abundance

# make base plot
plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "log(herbaceous flower abundance) (scaled)",
     ylab = "Pr(persistence)")

# add response curves for a range of species specialization metrics
for(i in 1:length(seq)){
  curve(ilogit(
    # phi0 +
    fit_summary$summary[19,1] + 
      # delta1_phi0*d[i] +
      fit_summary$summary[20,1]*seq[i] + 
      # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
      ((fit_summary$summary[22,1] + fit_summary$summary[23,1]*seq[i]) * x) +
      # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
      ((fit_summary$summary[25,1] + fit_summary$summary[26,1]*seq[i]) * 0)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # phi0 +
  fit_summary$summary[19,1] + 
    # delta1_phi0*d[i] +
    fit_summary$summary[20,1]*0 + 
    # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[22,1] + fit_summary$summary[23,1]*0) * x) +
    # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[25,1] + fit_summary$summary[26,1]*0) * 0)
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
    # phi0 +
    fit_summary$summary[19,1] + 
      # delta1_phi0*d[i] +
      fit_summary$summary[20,1]*seq[i] + 
      # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
      ((fit_summary$summary[22,1] + fit_summary$summary[23,1]*seq[i]) * 0) +
      # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
      ((fit_summary$summary[25,1] + fit_summary$summary[26,1]*seq[i]) * x)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # phi0 +
  fit_summary$summary[19,1] + 
    # delta1_phi0*d[i] +
    fit_summary$summary[20,1]*0 + 
    # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[22,1] + fit_summary$summary[23,1]*0) * 0) +
    # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[25,1] + fit_summary$summary[26,1]*0) * x)
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
    # phi0 +
    fit_summary$summary[1,1] + 
      # delta1_phi0*d[i] +
      fit_summary$summary[2,1]*seq[i] + 
      # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
      ((fit_summary$summary[4,1] + fit_summary$summary[5,1]*seq[i]) * x) +
      # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
      ((fit_summary$summary[7,1] + fit_summary$summary[8,1]*seq[i]) * 0)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # phi0 +
  fit_summary$summary[1,1] + 
    # delta1_phi0*d[i] +
    fit_summary$summary[2,1]*0 + 
    # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[4,1] + fit_summary$summary[5,1]*0) * x) +
    # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[7,1] + fit_summary$summary[8,1]*0) * 0)
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
    # phi0 +
    fit_summary$summary[1,1] + 
      # delta1_phi0*d[i] +
      fit_summary$summary[2,1]*seq[i] + 
      # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
      ((fit_summary$summary[4,1] + fit_summary$summary[5,1]*seq[i]) * 0) +
      # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
      ((fit_summary$summary[7,1] + fit_summary$summary[8,1]*seq[i]) * x)
  ),
  add=TRUE, col = my_palette[i], lwd = 3)
}

# add mean response curve (average d)
curve(ilogit(
  # phi0 +
  fit_summary$summary[1,1] + 
    # delta1_phi0*d[i] +
    fit_summary$summary[2,1]*0 + 
    # (delta0_phi_herbaceous + delta1_phi_herbaceous*d[i]) * herbaceous_flowers_scaled[x]
    ((fit_summary$summary[4,1] + fit_summary$summary[5,1]*0) * 0) +
    # hold at constant meany woody # (delta0_phi_woody + delta1_phi_woody*d[i]) * woody_flowers_scaled[] <- 0 
    ((fit_summary$summary[7,1] + fit_summary$summary[8,1]*0) * x)
),
add=TRUE, col = "black", lwd = 3, lty = "dashed")

