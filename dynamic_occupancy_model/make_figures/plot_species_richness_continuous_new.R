library(rstan)
library(viridis)
library(gridExtra)

source("./dynamic_occupancy_model/run_model/prep_data.R")
min_unique_detections = 1 # >=
my_data <- process_raw_data(min_unique_detections)

stan_out <- readRDS("./dynamic_occupancy_model/model_outputs/stan_out.rds")

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

original_herb_abundance <- my_data$original_herb_abundance
original_woody_abundance <- my_data$original_woody_abundance

## --------------------------------------------------
### Pull out model output

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## Use colour palette that is consistent with palette used for other figures

n_bins = 4 # 
n_years = 3
n_years_minus1 = n_years - 1
my_palette <- palette(viridis(n = n_bins, option = "viridis"))

## --------------------------------------------------
## get prediction range

n_draws = 2000 # number of samples from the posteriors

pred_length = 2

original_herb <- as.numeric(0:1) # reasonable prediction range 
herb_pred <- original_herb


#original_woody <- seq(min(original_woody_abundance), 
#                      max(original_woody_abundance), 
#                      length.out = pred_length) # reasonable prediction range
#woody_pred <- (original_woody - mean(original_woody_abundance)) / sd(original_woody_abundance) # unscale the dates

original_woody <- as.numeric(c(0,1)) # reasonable prediction range 
woody_pred <- original_woody

## --------------------------------------------------
## versus herbaceous flower abundance

psi1_expected <- array(data = NA, dim=c(n_species, pred_length))
gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years))

# expected values for a range of species and sites with range of herbaceous flowers
# BUT with average woody flower
random_draws_from_posterior = sample(1:nrow(list_of_draws), n_draws) 

richness = array(data = NA, dim=c( pred_length, n_years, n_draws))
psi <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
turnover <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1, n_draws))

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 2:n_years){
          
        psi1_expected[i,j] =
            ilogit(#psi1_0 +
              list_of_draws[rand,1] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + 4 + i
                list_of_draws[rand,(145+(i-1))] +
                #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
                list_of_draws[rand,3] * herb_pred[j] + 
                #psi1_specialization * d[i] +
                list_of_draws[rand,5] * d[i] + 
                #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
                list_of_draws[rand,6] * d[i] * herb_pred[j] 
            )
        
        gamma_expected[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,8] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(251+(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,10] * herb_pred[j] + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,12] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,13] * d[i] * herb_pred[j] 
          )
        
        phi_expected[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,15] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(357+(i-1))] +
              #psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
              list_of_draws[rand,17] * herb_pred[j] + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,19] * d[i] + 
              #(psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
              list_of_draws[rand,20] * d[i] * herb_pred[j]
          )
        
        if(k == 2){
          # calculate occurrence based on transition from first year
          psi[i,j,k-1] = psi1_expected[i,j] * phi_expected[i,j,k-1] + 
            (1 - psi1_expected[i,j]) * gamma_expected[i,j,k-1]
          
          turnover[i,j,k-1,draw] <- 
            ((1-psi1_expected[i,j]) * gamma_expected[i,j,k-1]) / 
            psi[i,j,1]
          
        } else{
          # calculate occurrence based on transition from first year
          psi[i,j,k-1] = psi[i,j,k-2] * phi_expected[i,j,k-1] + 
            (1 - psi[i,j,k-2]) * gamma_expected[i,j,k-1]
          
          turnover[i,j,k-1,draw] <- 
            ((1-psi[i,j,k-2]) * gamma_expected[i,j,k-1]) / 
            psi[i,j,2]
        }
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

df1 <- as.data.frame(cbind((rep(1, pred_length)), herb_pred,
                          mean[,1],
                          lower_50[,1], upper_50[,1],
                          lower_95[,1], upper_95[,1],
                          original_herb)) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df2 <- as.data.frame(cbind((rep(2, pred_length)), herb_pred,
                           mean[,2],
                           lower_50[,2], upper_50[,2],
                           lower_95[,2], upper_95[,2],
                           original_herb))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df3 <- as.data.frame(cbind((rep(3, pred_length)), herb_pred, 
                           mean[,3],
                           lower_50[,3], upper_50[,3],
                           lower_95[,3], upper_95[,3],
                           original_herb))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df <- rbind(df1, df2, df3) %>%
  mutate(year = as.factor(year),
         herb_pred = as.factor(herb_pred))


## --------------------------------------------------
## Draw species richness plot

(p  <- ggplot(data = df, aes(year, mean, colour=as.factor(original_herb))) +
   scale_x_discrete(name="Year", breaks = c(1, 2, 3),
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
   scale_colour_manual(name = "",
                     labels=c("control",
                              "herb. enhancement"),
                     values=my_palette[3:4]) +
   ylim(c(15, 110)) +
   theme_bw() +
   ylab("Expected species richness") +
   theme(legend.text=element_text(size=16),
         legend.position = c(0.75, 0.125),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

## --------------------------------------------------
## turnover summary

mean_t = matrix(nrow=pred_length, ncol=n_years_minus1)
lower_50_t = matrix(nrow=pred_length, ncol=n_years_minus1)
upper_50_t = matrix(nrow=pred_length, ncol=n_years_minus1)
lower_95_t = matrix(nrow=pred_length, ncol=n_years_minus1)
upper_95_t = matrix(nrow=pred_length, ncol=n_years_minus1)

for(j in 1:pred_length){
  for(k in 1:n_years_minus1){
    quants_t = as.vector(quantile(turnover[,j,k,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean_t[j,k] = quants_t[3]
    lower_50_t[j,k] = quants_t[2]
    upper_50_t[j,k] = quants_t[4]
    lower_95_t[j,k] = quants_t[1]
    upper_95_t[j,k] = quants_t[5]
  }
}

df1_t <- as.data.frame(cbind((rep(1, pred_length)), herb_pred,
                           mean_t[,1],
                           lower_50_t[,1], upper_50_t[,1],
                           lower_95_t[,1], upper_95_t[,1],
                           original_herb)) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df2_t <- as.data.frame(cbind((rep(2, pred_length)), herb_pred,
                           mean_t[,1],
                           lower_50_t[,1], upper_50_t[,1],
                           lower_95_t[,1], upper_95_t[,1],
                           original_herb))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df_t <- rbind(df1_t, df2_t) %>%
  mutate(year = as.factor(year),
         herb_pred = as.factor(herb_pred))


## --------------------------------------------------
## Draw turnover plot

(r  <- ggplot(data = df_t, aes(year, mean, colour=as.factor(original_herb))) +
   scale_x_discrete(name="Year", breaks = c(1, 2),
                    labels=c("2021/2022", "2022/2023"
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
   scale_colour_manual(name = "",
                       labels=c("control",
                                "herb. enhancement"),
                       values=my_palette[3:4]) +
   ylim(c(0, 1)) +
   theme_bw() +
   ylab("Average turnover rate") +
   theme(
         legend.position = "none",
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


psi1_expected_woody <- array(data = NA, dim=c(n_species, pred_length))
gamma_expected_woody <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
phi_expected_woody <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed_woody <- array(data = NA, dim=c(n_species, pred_length, n_years))

# expected values for a range of species and sites with range of herbaceous flowers
# BUT with average woody flower
random_draws_from_posterior = sample(1:nrow(list_of_draws), n_draws) 

richness_woody = array(data = NA, dim=c( pred_length, n_years, n_draws))
psi_woody <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
turnover_woody <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1, n_draws))

# this is set up to predict standardized for site having an herbaceous enhancement
for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 2:n_years){
        
        psi1_expected_woody[i,j] =
          ilogit(#psi1_0 +
            list_of_draws[rand,1] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(145+(i-1))] +
              # for an herb enahcnement site
              list_of_draws[rand,3] * 1 + 
              #psi1_specialization * d[i] +
              list_of_draws[rand,5] * d[i] + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,4] * woody_pred[j] + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,7] * d[i] * woody_pred[j]
          )
        
        gamma_expected_woody[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,8] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(251+(i-1))] +
              # for an herb enahcnement site
              list_of_draws[rand,10] * 1 +
              #psi1_specialization * d[i] +
              list_of_draws[rand,12] * d[i] + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,11] * woody_pred[j] + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,14] * d[i] * woody_pred[j]
          )
        
        phi_expected_woody[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,15] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + 4 + i
              list_of_draws[rand,(357+(i-1))] +
              # for an herb enahcnement site
              list_of_draws[rand,17] * 1 +
              #psi1_specialization * d[i] +
              list_of_draws[rand,19] * d[i] + 
              #psi1_woody_flowers * woody_flowers_scaled[j,1] +
              list_of_draws[rand,18] * woody_pred[j] + 
              #(psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
              list_of_draws[rand,21] * d[i] * woody_pred[j]
          )
        
        if(k == 2){
          # calculate occurrence based on transition from first year
          psi_woody[i,j,k-1] = psi1_expected_woody[i,j] * phi_expected_woody[i,j,k-1] + 
            (1 - psi1_expected_woody[i,j]) * gamma_expected_woody[i,j,k-1]
          
          turnover_woody[i,j,k-1,draw] <- 
            ((1-psi1_expected_woody[i,j]) * gamma_expected_woody[i,j,k-1]) / 
            psi_woody[i,j,1]
          
        } else{
          # calculate occurrence based on transition from first year
          psi_woody[i,j,k-1] = psi_woody[i,j,k-2] * phi_expected_woody[i,j,k-1] + 
            (1 - psi_woody[i,j,k-2]) * gamma_expected_woody[i,j,k-1]
          
          turnover_woody[i,j,k-1,draw] <- 
            ((1-psi_woody[i,j,k-2]) * gamma_expected_woody[i,j,k-1]) / 
            psi_woody[i,j,2]
        }
      } 
    }
  }
  
  # simmed occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 1:n_years){
        if(k < 2){
          occurrence_simmed_woody[i,j,1] <- rbinom(1, 1, prob = psi1_expected_woody[i,j])
        } else{
          occurrence_simmed_woody[i,j,k] = occurrence_simmed_woody[i,j,k-1] * phi_expected_woody[i,j,k-1] + 
            (1 - occurrence_simmed_woody[i,j,k-1]) * gamma_expected_woody[i,j,k-1]
        }
      }
    }
  }
  
  for(j in 1:pred_length){
    for(k in 1:n_years){
      richness_woody[j,k,draw] <- sum(occurrence_simmed_woody[1:n_species,j,k])
    }
  }
  
} 

mean_woody = matrix(nrow=pred_length, ncol=n_years)
lower_50_woody = matrix(nrow=pred_length, ncol=n_years)
upper_50_woody = matrix(nrow=pred_length, ncol=n_years)
lower_95_woody = matrix(nrow=pred_length, ncol=n_years)
upper_95_woody = matrix(nrow=pred_length, ncol=n_years)

for(j in 1:pred_length){
  for(k in 1:n_years){
    quants_woody = as.vector(quantile(richness_woody[j,k,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean_woody[j,k] = quants_woody[3]
    lower_50_woody[j,k] = quants_woody[2]
    upper_50_woody[j,k] = quants_woody[4]
    lower_95_woody[j,k] = quants_woody[1]
    upper_95_woody[j,k] = quants_woody[5]
  }
}

df1_woody <- as.data.frame(cbind((rep(1, pred_length)), woody_pred,
                           mean_woody[,1],
                           lower_50_woody[,1], upper_50_woody[,1],
                           lower_95_woody[,1], upper_95_woody[,1],
                           original_woody)) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 


df2_woody <- as.data.frame(cbind((rep(2, pred_length)), woody_pred,
                                 mean_woody[,1],
                                 lower_50_woody[,1], upper_50_woody[,1],
                                 lower_95_woody[,1], upper_95_woody[,1],
                                 original_woody))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df3_woody <- as.data.frame(cbind((rep(3, pred_length)), woody_pred,
                                 mean_woody[,1],
                                 lower_50_woody[,1], upper_50_woody[,1],
                                 lower_95_woody[,1], upper_95_woody[,1],
                                 original_woody))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 


df_woody <- rbind(df1_woody, df2_woody, df3_woody) %>%
  mutate(year = as.factor(year),
         woody_pred = as.factor(woody_pred))

## --------------------------------------------------
## Draw species richness plot

(q  <- ggplot(data = df_woody, aes(year, mean, colour=as.factor(original_woody))) +
   scale_x_discrete(name="Year", breaks = c(1, 2, 3),
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
   scale_colour_manual(name = "",
                       labels=c("avg. woody enhancement",
                                "above avg. \nwoody enhancement"),
                       values=my_palette[1:2]) +
   ylim(c(15, 110)) +
   theme_bw() +
   ylab("Expected species richness") +
   theme(legend.text=element_text(size=16),
         legend.position = c(0.7, 0.125),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

## --------------------------------------------------
## turnover summary

mean_t_woody = matrix(nrow=pred_length, ncol=n_years_minus1)
lower_50_t_woody = matrix(nrow=pred_length, ncol=n_years_minus1)
upper_50_t_woody = matrix(nrow=pred_length, ncol=n_years_minus1)
lower_95_t_woody = matrix(nrow=pred_length, ncol=n_years_minus1)
upper_95_t_woody = matrix(nrow=pred_length, ncol=n_years_minus1)

for(j in 1:pred_length){
  for(k in 1:n_years_minus1){
    quants_t_woody = as.vector(quantile(turnover_woody[,j,k,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean_t_woody[j,k] = quants_t_woody[3]
    lower_50_t_woody[j,k] = quants_t_woody[2]
    upper_50_t_woody[j,k] = quants_t_woody[4]
    lower_95_t_woody[j,k] = quants_t_woody[1]
    upper_95_t_woody[j,k] = quants_t_woody[5]
  }
}

df1_t_woody <- as.data.frame(cbind((rep(1, pred_length)), woody_pred,
                             mean_t_woody[,1],
                             lower_50_t_woody[,1], upper_50_t_woody[,1],
                             lower_95_t_woody[,1], upper_95_t_woody[,1],
                             original_woody)) %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df2_t_woody <- as.data.frame(cbind((rep(2, pred_length)), woody_pred,
                             mean_t_woody[,1],
                             lower_50_t_woody[,1], upper_50_t_woody[,1],
                             lower_95_t_woody[,1], upper_95_t_woody[,1],
                             original_woody))  %>%
  rename("year" = "V1",
         "mean" = "V3",
         "lower_50" = "V4",
         "upper_50" = "V5",
         "lower_95" = "V6",
         "upper_95" = "V7") 

df_t_woody <- rbind(df1_t_woody, df2_t_woody) %>%
  mutate(year = as.factor(year),
         woody_pred = as.factor(woody_pred))


## --------------------------------------------------
## Draw turnover plot

(s  <- ggplot(data = df_t_woody, aes(year, mean, colour=as.factor(original_woody))) +
   scale_x_discrete(name="Year", breaks = c(1, 2),
                    labels=c("2021/2022", "2022/2023"
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
   scale_colour_manual(name = "",
                       labels=c("avg. woody enhancement",
                                "above avg. \nwoody enhancement"),
                       values=my_palette[1:2]) +
   ylim(c(0, 1)) +
   theme_bw() +
   ylab("Average turnover rate") +
   theme(
     legend.position = "none",
     axis.text.x = element_text(size = 18),
     axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
     axis.title.x = element_text(size=20),
     axis.title.y = element_text(size = 20),
     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"))
)


#-------------------------------------------------------------------------------
grid.arrange(p, q, ncol=2)

#-------------------------------------------------------------------------------
grid.arrange(r, s, ncol=2)

#-------------------------------------------------------------------------------
grid.arrange(p, q, r, s, ncol=2)
