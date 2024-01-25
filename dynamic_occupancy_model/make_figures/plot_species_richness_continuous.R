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
my_palette_reduced <- rev(c(my_palette[2], my_palette[5], my_palette[7])) # purples
my_palette_reduced <- rev(c(my_palette[9], my_palette[12], my_palette[14])) # oranges

## --------------------------------------------------
## continuous habitat species richness plot

pred_seq <- seq(-2, 2, 0.1)
pred_length <- length(pred_seq)

psi1_expected <- array(data = NA, dim=c(n_species, pred_length))
gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years))

# expected values for a range of species and sites with range of herbaceous flowers
# BUT with average woody flower
n_draws = 100
random_draws_from_posterior = sample(1:nrow(list_of_draws), n_draws) 

richness = array(data = NA, dim=c( pred_length, n_years, n_draws))

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 2:n_years){
          
        psi1_expected[i,j] =
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
        
        gamma_expected[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,27] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(110+4*(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,28] * pred_seq[j] + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,30] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,31] * d[i] * pred_seq[j] + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,29] * 0 + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,32] * d[i] * 0 
          )
        
        phi_expected[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,33] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(111+4*(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,34] * pred_seq[j] + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,36] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,37] * d[i] * pred_seq[j] + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,35] * 0 + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,38] * d[i] * 0 
          )
      } 
    }
  }
  
  # simmed occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 1:n_years){
        if(k < 2){
          occurrence_simmed[i,j,1] <- rbinom(1, 1, prob = psi1_expected[i,j])
        } else{
          occurrence_simmed[i,j,k] = occurrence_simmed[i,j,k-1] * phi_expected[i,j,k-1] + 
            (1 - occurrence_simmed[i,j,k-1]) * gamma_expected[i,j,k-1]
        }
      }
    }
  }
  
  for(j in 1:pred_length){
    for(k in 1:n_years){
      richness[j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k])
    }
  }
    
} 
  
mean = matrix(nrow=pred_length, ncol=n_years)
lower_50 = matrix(nrow=pred_length, ncol=n_years)
upper_50 = matrix(nrow=pred_length, ncol=n_years)
lower_95 = matrix(nrow=pred_length, ncol=n_years)
upper_95 = matrix(nrow=pred_length, ncol=n_years)

for(j in 1:pred_length){
  for(k in 1:n_years){
    quants = as.vector(quantile(richness[j,k,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[j,k] = quants[3]
    lower_50[j,k] = quants[2]
    upper_50[j,k] = quants[4]
    lower_95[j,k] = quants[1]
    upper_95[j,k] = quants[5]
  }
}

df1 <- as.data.frame(cbind(as.factor(rep(1, pred_length)), pred_seq,
                          mean[,1],
                          lower_50[,1], upper_50[,1],
                          lower_95[,1], upper_95[,1])) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df2 <- as.data.frame(cbind(as.factor(rep(2, pred_length)), pred_seq,
                           mean[,2],
                           lower_50[,2], upper_50[,2],
                           lower_95[,2], upper_95[,2]))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df3 <- as.data.frame(cbind(as.factor(rep(3, pred_length)), pred_seq,
                           mean[,3],
                           lower_50[,3], upper_50[,3],
                           lower_95[,3], upper_95[,3]))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df <- rbind(df1, df2, df3) 


## --------------------------------------------------
## Draw species richness plot

(p <- ggplot(data = df, aes(pred_seq, mean)) +
  geom_line(data=df1, size=2) +
  geom_ribbon(data=df1, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[1], alpha=0.8) + 
  geom_ribbon(data=df1, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[1], alpha=0.1) + 
  geom_line(data=df2, size=2) +
  geom_ribbon(data=df2, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[2], alpha=0.8) + 
  geom_ribbon(data=df2, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[2], alpha=0.1) + 
  geom_line(data=df3, size=2) +
  geom_ribbon(data=df3, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[3], alpha=0.8) + 
  geom_ribbon(data=df3, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[3], alpha=0.1) +  
  ylim(c(15, 85)) +
  theme_bw() +
  ylab("predicted species richness") +
  xlab("herbaceous flower abundance") +
  #scale_color_manual(values=my_palette_reduced) +
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

psi1_expected <- array(data = NA, dim=c(n_species, pred_length))
gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years))

# expected values for a range of species and sites with range of herbaceous flowers
# BUT with average woody flower
n_draws = 100
random_draws_from_posterior = sample(1:nrow(list_of_draws), n_draws) 

richness = array(data = NA, dim=c( pred_length, n_years, n_draws))

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 2:n_years){
        
        psi1_expected[i,j] =
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
        
        gamma_expected[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,27] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(110+4*(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,28] * 0 + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,30] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,31] * d[i] * 0 + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,29] * pred_seq[j] + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,32] * d[i] * pred_seq[j]
          )
        
        phi_expected[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,33] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(111+4*(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,34] * 0 + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,36] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,37] * d[i] * 0 + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,35] * pred_seq[j] + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,38] * d[i] * pred_seq[j]
          )
      } 
    }
  }
  
  # simmed occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 1:n_years){
        if(k < 2){
          occurrence_simmed[i,j,1] <- rbinom(1, 1, prob = psi1_expected[i,j])
        } else{
          occurrence_simmed[i,j,k] = occurrence_simmed[i,j,k-1] * phi_expected[i,j,k-1] + 
            (1 - occurrence_simmed[i,j,k-1]) * gamma_expected[i,j,k-1]
        }
      }
    }
  }
  
  for(j in 1:pred_length){
    for(k in 1:n_years){
      richness[j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k])
    }
  }
  
} 

mean = matrix(nrow=pred_length, ncol=n_years)
lower_50 = matrix(nrow=pred_length, ncol=n_years)
upper_50 = matrix(nrow=pred_length, ncol=n_years)
lower_95 = matrix(nrow=pred_length, ncol=n_years)
upper_95 = matrix(nrow=pred_length, ncol=n_years)

for(j in 1:pred_length){
  for(k in 1:n_years){
    quants = as.vector(quantile(richness[j,k,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[j,k] = quants[3]
    lower_50[j,k] = quants[2]
    upper_50[j,k] = quants[4]
    lower_95[j,k] = quants[1]
    upper_95[j,k] = quants[5]
  }
}

df1 <- as.data.frame(cbind(as.factor(rep(1, pred_length)), pred_seq,
                           mean[,1],
                           lower_50[,1], upper_50[,1],
                           lower_95[,1], upper_95[,1])) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df2 <- as.data.frame(cbind(as.factor(rep(2, pred_length)), pred_seq,
                           mean[,2],
                           lower_50[,2], upper_50[,2],
                           lower_95[,2], upper_95[,2]))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df3 <- as.data.frame(cbind(as.factor(rep(3, pred_length)), pred_seq,
                           mean[,3],
                           lower_50[,3], upper_50[,3],
                           lower_95[,3], upper_95[,3]))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df <- rbind(df1, df2, df3) 

## --------------------------------------------------
## Draw species richness plot
# use a new palette selection



(q  <- ggplot(data = df, aes(pred_seq, mean)) +
    geom_line(data=df1, size=2) +
    geom_ribbon(data=df1, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[1], alpha=0.8) + 
    geom_ribbon(data=df1, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[1], alpha=0.1) + 
    geom_line(data=df2, size=2) +
    geom_ribbon(data=df2, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[2], alpha=0.8) + 
    geom_ribbon(data=df2, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[2], alpha=0.1) + 
    geom_line(data=df3, size=2) +
    geom_ribbon(data=df3, aes(ymin=lower_50, ymax=upper_50), fill = my_palette_reduced[3], alpha=0.8) + 
    geom_ribbon(data=df3, aes(ymin=lower_95, ymax=upper_95), fill = my_palette_reduced[3], alpha=0.1) +  
    ylim(c(15, 85)) +
    theme_bw() +
    ylab("predicted species richness") +
    xlab("woody flower abundance") +
    #scale_color_manual(values=my_palette_reduced) +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)
