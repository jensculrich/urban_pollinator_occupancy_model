library(rstan)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 2 # >=
my_data <- process_raw_data(min_unique_detections)

stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out3_2ormore.rds")

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

list_of_draws <- as.data.frame(stan_out)
colnames(list_of_draws)

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
## continuous habitat species richness plot

pred_seq <- seq(-2, 2, 0.1)
pred_length <- length(pred_seq)

occurrence1_expected <- matrix(data = NA, nrow = n_species, ncol = pred_length)
occurrence1_simmed <- matrix(data = NA, nrow = n_species, ncol = pred_length)

# expected values for a range of species and sites with range of herbaceous flowers
# BUT with average woody flower
n_draws = 500
random_draws_from_posterior = sample(1:nrow(list_of_draws), n_draws) 

richness = matrix(nrow = pred_length, ncol = n_draws)

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  for(i in 1:n_species){
    for(j in 1:pred_length){
      occurrence1_expected[i,j] =
        ilogit(#psi1_0 +
          list_of_draws[rand,21] + 
            #species_effects[species[i],1] + // a species specific intercept
            # start at first row of species effects
            # then each next species will be + 4 + i
            list_of_draws[rand,(109+4*(i-1))] +
            #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
            list_of_draws[rand,22] * pred_seq[j] + 
            #psi1_specialization * d[i] +
            list_of_draws[rand,24] * d[i] + 
            #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
            list_of_draws[rand,25] * d[i] * pred_seq[j] + 
            #psi1_woody_flowers * woody_flowers_scaled[j,1] +
            list_of_draws[rand,23] * 0 + 
            #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
            list_of_draws[rand,26] * d[i] * 0 
        )
    }
  }
  
  for(i in 1:n_species){
    for(j in 1:pred_length){
      occurrence1_simmed[i,j] <- rbinom(1, 1, prob = occurrence1_expected[i,j])
    }
  }
  
  for(j in 1:pred_length){
    richness[j,draw] <- sum(occurrence1_simmed[1:n_species,j])
  }
}

mean = vector(length=pred_length)
lower_50 = vector(length=pred_length)
upper_50 = vector(length=pred_length)
lower_95 = vector(length=pred_length)
upper_95 = vector(length=pred_length)

for(i in 1:pred_length){
  quants = as.vector(quantile(richness[i,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
  
  mean[i] = quants[3]
  lower_50[i] = quants[2]
  upper_50[i] = quants[4]
  lower_95[i] = quants[1]
  upper_95[i] = quants[5]
  
}

df <- as.data.frame(cbind(pred_seq,
                          mean,
                          lower_50, upper_50,
                          lower_95, upper_95))

## --------------------------------------------------
## Draw species richness plot

(p <- ggplot(data = df, aes(pred_seq, mean)) +
  geom_line(size=2) +
  geom_ribbon(aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[1], alpha=0.8) + 
  geom_ribbon(aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[1], alpha=0.1) + 
  ylim(c(20, 85)) +
  theme_bw() +
  ylab("predicted species richness (year 1)") +
  xlab("herbaceous flower abundance") +
  scale_color_manual(values=my_palette_reduced) +
  theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)


##------------------------------------------------------------------------------
## Now do it with woody plant richness

## --------------------------------------------------
## continuous habitat species richness plot

pred_seq <- seq(-2, 2, 0.1)
pred_length <- length(pred_seq)

occurrence1_expected2 <- matrix(data = NA, nrow = n_species, ncol = pred_length)
occurrence1_simmed2 <- matrix(data = NA, nrow = n_species, ncol = pred_length)

# expected values for a range of species and sites with range of woody flowers
# BUT with average woody flower
n_draws = 500
random_draws_from_posterior2 = sample(1:nrow(list_of_draws), n_draws) 

richness2 = matrix(nrow = pred_length, ncol = n_draws)

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior2[draw]
  
  for(i in 1:n_species){
    for(j in 1:pred_length){
      occurrence1_expected2[i,j] =
        ilogit(#psi1_0 +
          list_of_draws[rand,21] + 
            #species_effects[species[i],1] + // a species specific intercept
            # start at first row of species effects
            # then each next species will be + 4 + i
            list_of_draws[rand,(109+4*(i-1))] +
            #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
            list_of_draws[rand,22] * 0 + 
            #psi1_specialization * d[i] +
            list_of_draws[rand,24] * d[i] + 
            #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
            list_of_draws[rand,25] * d[i] * 0 + 
            #psi1_woody_flowers * woody_flowers_scaled[j,1] +
            list_of_draws[rand,23] * pred_seq[j] + 
            #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
            list_of_draws[rand,26] * d[i] * pred_seq[j]
        )
    }
  }
  
  for(i in 1:n_species){
    for(j in 1:pred_length){
      occurrence1_simmed2[i,j] <- rbinom(1, 1, prob = occurrence1_expected2[i,j])
    }
  }
  
  for(j in 1:pred_length){
    richness2[j,draw] <- sum(occurrence1_simmed2[1:n_species,j])
  }
}

mean2 = vector(length=pred_length)
lower_50_2 = vector(length=pred_length)
upper_50_2 = vector(length=pred_length)
lower_95_2 = vector(length=pred_length)
upper_95_2 = vector(length=pred_length)

for(i in 1:pred_length){
  quants = as.vector(quantile(richness2[i,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
  
  mean2[i] = quants[3]
  lower_50_2[i] = quants[2]
  upper_50_2[i] = quants[4]
  lower_95_2[i] = quants[1]
  upper_95_2[i] = quants[5]
  
}

df2 <- as.data.frame(cbind(pred_seq,
                          mean2,
                          lower_50_2, upper_50_2,
                          lower_95_2, upper_95_2))

## --------------------------------------------------
## Draw species richness plot

(p <- ggplot(data = df2, aes(pred_seq, mean2)) +
   geom_line(size=2) +
   geom_ribbon(aes(ymin=lower_50_2, ymax=upper_50_2), fill = my_palette_reduced[2], alpha=0.8) + 
   geom_ribbon(aes(ymin=lower_95_2, ymax=upper_95_2), fill = my_palette_reduced[2], alpha=0.1) + 
   ylim(c(20, 85)) +
   theme_bw() +
   ylab("predicted species richness (year 1)") +
   xlab("woody flower abundance") +
   scale_color_manual(values=my_palette_reduced) +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)
