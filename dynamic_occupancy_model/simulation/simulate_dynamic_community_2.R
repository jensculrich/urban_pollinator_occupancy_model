library(sn) # for multivariate skewed normal distr data simulation


## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_species <- 100 # number of species
n_sites <- 20 # number of sites (must be an even number for simulation code)
n_years <- 4 # number of years
n_years_minus1 <- n_years - 1
n_visits <- 6 # number of surveys per year

mean_enhancement_herb_flowers <- 0.5
sigma_herb_flowers_means <- 0.5
sigma_herbaceous_flowers_visit <- 1
sigma_woody_flowers_means <- 1
sigma_woody_flowers_visit <- 1

# set parameter values
psi1_0 <- 0 # prob of initial occupancy
sigma_psi1_species <- 2 # prob of initial occupancy
psi1_herbaceous_flowers <- 0.6
psi1_woody_flowers <- 0.4
psi1_specialization <- -1
psi1_interaction_1 <- 0
psi1_interaction_2 <- 0.5

gamma0 <- -0.25 # prob of initial occupancy
sigma_gamma_species <- 1 # prob of initial occupancy
gamma_herbaceous_flowers <- 0
gamma_woody_flowers <- 0
gamma_specialization <- -1.25
gamma_interaction_1 <- 0
gamma_interaction_2 <- 1
gamma_min = 0 # set these to zero if you want no year heterogeneity
gamma_max = 0 

phi0 <- 2 # prob of initial occupancy
sigma_phi_species <- 1 # prob of initial occupancy
phi_herbaceous_flowers <- 0
phi_woody_flowers <- 0
phi_specialization <- -0.5
phi_interaction_1 <- 0.5
phi_interaction_2 <- 0.5
phi_min = 0 # set these to zero if you want no year heterogeneity
phi_max = 0 

p0 <- 0 # probability of detection (logit scaled)
sigma_p_species <- 2 # species-specific variation
p_degree <- 0.75
p_flower_abundance_any <- 0.5 # increase in detection rate moving from one habitat type to the other (logit scaled)
mu_p_species_date <- 0
sigma_p_species_date <- 1
mu_p_species_date_sq <- -0.5  
sigma_p_species_date_sq <- 1
p_overdispersion_sigma <- 0.5

mean_survey_date <- 180
sigma_survey_date <- 40

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- FALSE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
    n_species, n_sites, n_years, n_years_minus1, n_visits,
    psi1_0,
    sigma_psi1_species,
    psi1_herbaceous_flowers,
    psi1_woody_flowers,
    psi1_specialization,
    psi1_interaction_1,
    psi1_interaction_2,
    
    gamma0,
    sigma_gamma_species,
    gamma_herbaceous_flowers,
    gamma_woody_flowers,
    gamma_specialization,
    gamma_interaction_1,
    gamma_interaction_2,
    
    phi0,
    sigma_phi_species,
    phi_herbaceous_flowers,
    phi_woody_flowers,
    phi_specialization,
    phi_interaction_1,
    phi_interaction_2,
    
    p0, # probability of detection (logit scaled)
    sigma_p_species, # species-specific variation
    p_degree,
    p_flower_abundance_any, # increase in detection rate moving from one habitat type to the other (logit scaled)
    mu_p_species_date,
    sigma_p_species_date,
    mu_p_species_date_sq,  
    sigma_p_species_date_sq,
    p_overdispersion_sigma,
    mean_survey_date,
    sigma_survey_date,
    create_missing_data,
    prob_missing
){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_species, n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_species, n_sites, n_years, n_visits)) # observed data
  
  gamma_year <- vector(length = n_years_minus1)
  for(k in 1:length(gamma_year)){
    gamma_year[k] <- runif(1, min=gamma_min, max=gamma_max)
  }
  
  phi_year <- vector(length = n_years_minus1)
  for(k in 1:length(phi_year)){
    phi_year[k] <- runif(1, min=phi_min, max=phi_max)
  }
  
  # p <- p # probability of detection
  # equilibrium occupancy rate (for the average species)
  #(psi_eq_habitat0 <- ilogit(gamma0) / (ilogit(gamma0)+(1-ilogit(phi0)))) # equilibrium occupancy rate 
  #(psi_eq_habitat1 <- ilogit(gamma0 + mu_gamma_habitat) / (ilogit(gamma0 + mu_gamma_habitat)+(1-ilogit(phi0 + mu_phi_habitat)))) # equilibrium occupancy rate
  
  ## --------------------------------------------------
  ## Create covariate data
  
  ## --------------------------------------------------
  ## species specialization
  
  #d <- rbeta(n_species, 0.75, 2)
  #d_scaled <- center_scale(d)
  
  # in my real data most species tend to have pretty low specialization (skewed left)
  # the degree and d' are possibly positively correlated
  
  xi <- c(0, 0)
  Omega <- diag(2)
  Omega[2,1] <- Omega[1,2] <- 0.5
  alpha <- c(3,3)
  correlated_data <- sn::rmsn(n_species, xi, Omega, alpha)
  
  d <- correlated_data[,1]
  degree <- correlated_data[,2]
  
  #define function to scale values between 0 and 1
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  
  #scale values in 'sales' column to be between 0 and 1
  d <- scale_values(d)
  degree <- scale_values(degree)
  
  #hist(d)
  #hist(degree)
  
  d_scaled <- center_scale(d)
  degree_scaled <- center_scale(degree)
  
  #hist(d_scaled)
  #hist(degree_scaled)
  
  ## --------------------------------------------------
  ## day of year
  # should have a value for each site*year*visit
  date <- array(NA, dim =c(n_sites, n_years, n_visits))
  mean_survey_date = mean_survey_date
  sigma_survey_date = sigma_survey_date
  
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      # create a vector of visit dates centered on the middle of the early summer
      date[j, k, 1:n_visits] <- sort(as.integer(rnorm(
        n_visits, mean = mean_survey_date, sd = sigma_survey_date))) 
    }
  }
  
  ## let's scale the calendar day of year by the mean date (z-score scaled)
  date_scaled <- center_scale(date) 
  
  ## --------------------------------------------------
  ## habitat type
  
  ## habitat type
  habitat_type <- rep(c(0,1), each = n_sites / 2)
  
  ## --------------------------------------------------
  ## herbaceous flower cover
  
  # flower abundance herb by survey
  
  # there is probably some connection between visits within each given site
  # so we will make some mean site values where some sites tend to have more flowers
  # for the herbs this is in part driven by the management
  herbaceous_flower_means <- vector(length=n_sites)
  for(j in 1:n_sites){
    herbaceous_flower_means[j] = rnorm(1, mean=mean_enhancement_herb_flowers*habitat_type[j], sd=sigma_herb_flowers_means) 
  }
  # View(as.data.frame(cbind(habitat_type, herbaceous_flower_means)))
  mean(herbaceous_flower_means[1:10])
  mean(herbaceous_flower_means[11:20])
  
  # now generate visit specific herbaceous flower counts
  herbaceous_flowers_scaled_by_survey <- array(data = NA, dim = c(n_sites, n_years, n_visits))
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      for(l in 1:n_visits){
        herbaceous_flowers_scaled_by_survey[j,k,l] = rnorm(1, herbaceous_flower_means[j], sigma_herbaceous_flowers_visit)
      }
    }
  }
  mean(herbaceous_flowers_scaled_by_survey[1:(n_sites/2),,])
  mean(herbaceous_flowers_scaled_by_survey[((n_sites/2)+1):n_sites,,])
  
  ## --------------------------------------------------
  ## woody flower cover
  
  # flower abundance woody by survey
  
  # there is probably some connection between visits within each given site
  # so we will make some mean site values where some sites tend to have more flowers
  # for the woody this is unrelated to the habitat management and pretty random
  woody_flower_means <- vector(length=n_sites)
  for(j in 1:n_sites){
    woody_flower_means[j] = rnorm(1, mean=0, sd=sigma_woody_flowers_means) 
  }
  # View(as.data.frame(cbind(habitat_type, herbaceous_flower_means)))
  mean(woody_flower_means[1:10])
  mean(woody_flower_means[11:20])
  
  # now generate visit specific herbaceous flower counts
  woody_flowers_scaled_by_survey <- array(data = NA, dim = c(n_sites, n_years, n_visits))
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      for(l in 1:n_visits){
        woody_flowers_scaled_by_survey[j,k,l] = rnorm(1, woody_flower_means[j], sigma_woody_flowers_visit)
      }
    }
  }
  mean(woody_flowers_scaled_by_survey[1:10,,])
  mean(woody_flowers_scaled_by_survey[11:20,,])
  
  
  # flower abundance woody by survey
  # these two things should be correlated
  # then calculate average annual of each (for ecological processes)
  # and calculate average of both of together for each survey (for detection)
    
  flowers_any_by_survey <- array(data = NA, dim = c(n_sites, n_years, n_visits))
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      for(l in 1:n_visits){
        flowers_any_by_survey[j,k,l] = 
          (herbaceous_flowers_scaled_by_survey[j,k,l] + woody_flowers_scaled_by_survey[j,k,l]) / 2
          
      }
    }
  }
  mean(flowers_any_by_survey[1:10,,])
  mean(flowers_any_by_survey[11:20,,])
  
  ## --------------------------------------------------
  ## calculate mean flowers by site X year for each category
  
  # herbaceous flowers
  herbaceous_flowers_scaled <- matrix(NA, nrow=n_sites, ncol=n_years)
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      herbaceous_flowers_scaled[j,k] = mean(herbaceous_flowers_scaled_by_survey[j,k,])
    }
  }
  
  # woody flowers
  woody_flowers_scaled <- matrix(NA, nrow=n_sites, ncol=n_years)
  for(j in 1:n_sites) { # for each site
    for(k in 1:n_years) { # for each n_years
      woody_flowers_scaled[j,k] = mean(woody_flowers_scaled_by_survey[j,k,])
    }
  }
  
  ## --------------------------------------------------
  ## Create random effects
  
  ## species-specific random intercepts
  psi1_species <- rnorm(n=n_species, mean=psi1_0, sd=sigma_psi1_species)
  
  gamma_species <- rnorm(n=n_species, mean=gamma0, sd=sigma_gamma_species)

  phi_species <- rnorm(n=n_species, mean=phi0, sd=sigma_phi_species)

  p_species <- rnorm(n=n_species, mean=p0, sd=sigma_p_species)

  # p_site <- rnorm(n=n_sites, mean=0, sd=sigma_p_site)
  
  p_species_date <- rnorm(n=n_species, mean=mu_p_species_date, sd=sigma_p_species_date)
  
  p_species_date_sq <- rnorm(n=n_species, mean=mu_p_species_date_sq, sd=sigma_p_species_date_sq)
  
  ## --------------------------------------------------
  ## Create expected values
  
  # generate p with heterogeneity
  logit_p <- array(NA, dim = c(n_species, n_sites, n_years, n_visits)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_visits){
          
          logit_p[i,j,k,l] = 
            p_species[i] +
            p_degree * degree[i] +
            p_flower_abundance_any * flowers_any_by_survey[j,k,l] +
            p_species_date[i]*date_scaled[j, k, l] + # a spatiotemporally specific intercept
            p_species_date_sq[i]*(date_scaled[j, k, l])^2 # a spatiotemporally specific intercept
          
        }
      }
    }  
  }

  
  # generate ecological expected values with heterogeneity
  logit_psi1 <- array(NA, dim = c(n_species, n_sites)) 
  logit_gamma <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  logit_phi <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years_minus1){
        
        logit_psi1[i,j] = 
          psi1_species[i] +
          psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,k] +  
          psi1_specialization * d[i] +
          (psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,k]) +
          psi1_woody_flowers * woody_flowers_scaled[j,k] + 
          (psi1_interaction_2 * d[i] * woody_flowers_scaled[j,k])
        
        logit_gamma[i,j,k] = 
          gamma_species[i] +
          gamma_herbaceous_flowers * herbaceous_flowers_scaled[j,k] + 
          gamma_specialization * d[i] +
          (gamma_interaction_1 * d[i] * herbaceous_flowers_scaled[j,k]) +
          gamma_woody_flowers * woody_flowers_scaled[j,k] + 
          (gamma_interaction_2 * d[i] * woody_flowers_scaled[j,k])
        
        logit_phi[i,j,k] = 
          phi_species[i] +
          phi_herbaceous_flowers * herbaceous_flowers_scaled[j,k] + 
          phi_specialization * d[i] +
          (phi_interaction_1 * d[i] * herbaceous_flowers_scaled[j,k]) +
          phi_woody_flowers * woody_flowers_scaled[j,k] + 
          (phi_interaction_2 * d[i] * woody_flowers_scaled[j,k])
      
      }
    }
  }
        


  # generate initial presence/absence states
  for(i in 1:n_species){
    for(j in 1:n_sites){
      z[i,j,1] <- rbinom(n=1, size=1, prob=ilogit(logit_psi1[i,j])) 
    }
  }

  sum(z[,,1]) / n_sites # true occupancy proportion in year 1
  
  # generate presence/absence in subsequent years
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 2:n_years){
        
        # use z as a switch so we are estimating 
        exp_z <- z[i,j,k-1] * ilogit(logit_phi[i,j,k-1]) + # survival if z=1
          (1 - z[i,j,k-1]) * ilogit(logit_gamma[i,j,k-1]) # or colonization if z=0
        
        # and then assign z stochastically
        # some sites may transition if they are colonized or local extinction occurs
        # but might otherwise retain their state across years
        z[i,j,k] <- rbinom(n = 1, size = 1, prob = exp_z) 
        
      }
    }    
  }

  apply(z, 3, sum) / n_sites # true occupancy proportions in each year
  
  apply(z[,1:(n_sites/2),], 3, sum) / (n_sites / 2) # true occupancy proportions in each year
  apply(z[,(n_sites/2+1):n_sites,], 3, sum) / (n_sites / 2) # true occupancy proportions in each year
  
  # detection / non-detection data
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_visits){
          y[i,j,k,l] <- rbinom(n = 1, size = 1, prob = z[i,j,k]*ilogit(logit_p[i,j,k,l]))
        }
      }
    }
  }

  y ; str(y)
  
  sum((y/n_visits) / sum(z)) # proportion of times detection given presence
  
  # need to fix this if going to use
  if(create_missing_data == TRUE){
    # generate missing values: create simple version of unbalanced data set
    prob_missing <- prob_missing # constant NA probability
    y2 <- y # duplicate balanced dataset
    for(i in 1:n_sites){
      for(j in 1:n_visits){
        for(t in 1:n_years){
          turnNA <- rbinom(1,1,prob_missing)
          y2[i,t,j] <- ifelse(turnNA==1, NA, y2[i,t,j])
        }
      }
    }
    y2 ; str(y2)
  } else {
    y2 <- y # duplicate dataset with no missing data
  }
  
  # which species never occurred
  species_not_occurring <- length(which(apply(z, 1, sum) == 0))
  
  # which species were never detected (even though they occurred)
  species_not_observed <- length(which(apply(y2, 1, sum) == 0))
  
  table(nvisits <- apply(y2, c(1,3), function(x) sum(!is.na(x))))
  # _ sites with no visits, _ with 1 visit, and _ with 2 visits through 7 visits
  
  # compute true expected and realized occupancy (psi and psi_fs)
  #psi_habitat0 <- numeric(n_years) ; psi_habitat0[1] <- ilogit(psi1_0)
  #psi_habitat1 <- numeric(n_years) ; psi_habitat1[1] <- ilogit(psi1_0 + mu_psi1_habitat)
  
  # expected for the average species across the two habitat types
    #for(k in 2:n_years){ # compute true values of psi
    #  psi_habitat0[k] <- psi_habitat0[k-1] * ilogit(phi0) + (1 - psi_habitat0[k-1]) * ilogit(gamma0)
    #}
  
    #for(k in 2:n_years){ # compute true values of psi
   #   psi_habitat1[k] <- psi_habitat1[k-1] * ilogit(phi0 + mu_phi_habitat) + (1 - psi_habitat1[k-1]) * ilogit(gamma0 + mu_gamma_habitat)
   # }
  
  # actual occupancy rate given the simulated data (for the average species)
  # only works if 1 species # psi_fs_habitat0 <- colSums(z[1:n_species,1:(n_sites/2),]) / (n_sites/2) # psi finite sample
  #psi_fs_habitat0 <- apply(z[,1:(n_sites/2),], 3, function(x) mean(x))
  # only works if 1 species # psi_fs_habitat1 <- colSums(z[,(n_sites/2+1):n_sites,]) / (n_sites/2) # psi finite sample
  #psi_fs_habitat1 <- apply(z[,(n_sites/2+1):n_sites,], 3, function(x) mean(x))
  
  # compute observed occupancy proportion (for the average species)
  zobs <- apply(y2, c(1,2,3), function(x) max(x, na.rm = TRUE))
  zobs[zobs == "-Inf"] <- NA
  psi_obs_habitat0 <- apply(zobs[,1:(n_sites/2),], 3, sum, na.rm = TRUE) / apply(zobs, 3, function(x) sum(!is.na(x)))
  psi_obs_habitat1 <- apply(zobs[,(n_sites/2+1):n_sites,], 3, sum, na.rm = TRUE) / apply(zobs, 3, function(x) sum(!is.na(x)))
  
  print(paste0("*** ", species_not_occurring, " / ", n_species, " species never occurred ***"))
  print(paste0("*** ", species_not_observed, " / ", n_species, " species never detected ***"))
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    #psi_habitat0 = psi_habitat0,
    #psi_habitat1 = psi_habitat1,
    #psi_fs_habitat0 = psi_fs_habitat0,
    #psi_fs_habitat1 = psi_fs_habitat1,
    #psi_obs_habitat0 = psi_obs_habitat0,
    #psi_obs_habitat1 = psi_obs_habitat1,
    #psi_eq_habitat0 = psi_eq_habitat0,
    #psi_eq_habitat1 = psi_eq_habitat1,
    V = y2, # return detection data after potentially introducing NAs,
    habitat_type = habitat_type,
    herbaceous_flowers_scaled = herbaceous_flowers_scaled,
    woody_flowers_scaled = woody_flowers_scaled,
    flowers_any_by_survey = flowers_any_by_survey,
    
    date_scaled = date_scaled,
    gamma_year = gamma_year,
    phi_year = phi_year,
    d = d_scaled,
    degree = degree_scaled,
    flowers_any_by_survey
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

set.seed(1)
my_simulated_data <- simulate_data(  
  n_species, n_sites, n_years, n_years_minus1, n_visits,
  psi1_0,
  sigma_psi1_species,
  psi1_herbaceous_flowers,
  psi1_woody_flowers,
  psi1_specialization,
  psi1_interaction_1,
  psi1_interaction_2,
  
  gamma0,
  sigma_gamma_species,
  gamma_herbaceous_flowers,
  gamma_woody_flowers,
  gamma_specialization,
  gamma_interaction_1,
  gamma_interaction_2,
  
  phi0,
  sigma_phi_species,
  phi_herbaceous_flowers,
  phi_woody_flowers,
  phi_specialization,
  phi_interaction_1,
  phi_interaction_2,
  
  p0, # probability of detection (logit scaled)
  sigma_p_species, # species-specific variation
  p_degree,
  p_flower_abundance_any, # increase in detection rate moving from one habitat type to the other (logit scaled)
  mu_p_species_date,
  sigma_p_species_date,
  mu_p_species_date_sq,  
  sigma_p_species_date_sq,
  p_overdispersion_sigma,
  mean_survey_date,
  sigma_survey_date,
  create_missing_data,
  prob_missing
  )

V <- my_simulated_data$V
habitat_type <- my_simulated_data$habitat_type
#psi_habitat0 <- my_simulated_data$psi_habitat0
#psi_habitat1 <- my_simulated_data$psi_habitat1
#psi_fs_habitat0 <- my_simulated_data$psi_fs_habitat0
#psi_fs_habitat1 <- my_simulated_data$psi_fs_habitat1
#psi_obs_habitat0 <- my_simulated_data$psi_obs_habitat0
#psi_obs_habitat1 <- my_simulated_data$psi_obs_habitat1
#psi_eq_habitat0 <- my_simulated_data$psi_eq_habitat0
#psi_eq_habitat1 <- my_simulated_data$psi_eq_habitat1
date_scaled <- my_simulated_data$date_scaled
d <- my_simulated_data$d
degree <- my_simulated_data$degree
herbaceous_flowers_scaled <- my_simulated_data$herbaceous_flowers_scaled
woody_flowers_scaled <- my_simulated_data$woody_flowers_scaled
flowers_any_by_survey <- my_simulated_data$flowers_any_by_survey
species <- seq(1, n_species, by=1)
sites <- seq(1, n_sites, by=1)
years <- seq(1, n_years_minus1, by=1)

gamma_year <- my_simulated_data$gamma_year
phi_year <- my_simulated_data$phi_year

## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "species", "sites", "years",
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_visits",
               "habitat_type", "date_scaled", "d", "degree", 
               "herbaceous_flowers_scaled", "woody_flowers_scaled", "flowers_any_by_survey") 


## Parameters monitored 
params <- c(
  "psi1_0", "sigma_psi1_species",
  "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
  "psi1_interaction_1", "psi1_interaction_2",
  
  "gamma0", "sigma_gamma_species",
  "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
  "gamma_interaction_1", "gamma_interaction_2",
  
  "phi0", "sigma_phi_species",
  "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
  "phi_interaction_1", "phi_interaction_2",
  
  "p0", "sigma_p_species", "p_degree", 
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq", "p_flower_abundance_any", 
  "species_richness", "avg_species_richness_control", "avg_species_richness_enhanced", "increase_richness_enhanced",
  #"turnover_control", "turnover_enhanced",
  #"psi_eq_habitat0", "psi_eq_habitat1",
  "T_rep", "T_obs", "P_species")

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains

# targets
parameter_values <-  c(
  psi1_0, sigma_psi1_species, mu_psi1_habitat, sigma_psi1_habitat,
  gamma0, sigma_gamma_species, mu_gamma_habitat, sigma_gamma_habitat, NA,
  phi0, sigma_phi_species, mu_phi_habitat, sigma_phi_habitat, NA,
  p0, sigma_p_species, sigma_p_site, p_habitat, 
  mu_p_species_date, sigma_p_species_date, mu_p_species_date_sq, sigma_p_species_date_sq,
  NA, NA, NA
)



## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0, 1),
       mu_psi1_habitat = runif(1, -1, 1),
       sigma_psi1_habitat = runif(1, 0, 1),
       gamma0 = runif(1, -1, 1),
       sigma_gamma_species = runif(1, 0, 1),
       mu_gamma_habitat = runif(1, -2, -1), # colonization rates are usually low
       sigma_gamma_habitat = runif(1, 0, 1),
       phi0 = runif(1, -1, 1),
       sigma_phi_species = runif(1, 0, 1),
       mu_phi_habitat = runif(1, 0, 1), # persistence rates are usually greater than 50%
       sigma_phi_habitat = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
       p_habitat = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

targets <- as.data.frame(cbind(params, parameter_values))


## --------------------------------------------------
### Run model

library(rstan)
stan_model <- "./dynamic_occupancy_model/models/dynocc_model_4.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     #init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

saveRDS(stan_out_sim, "./dynamic_occupancy_model/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./dynamic_occupancy_model/simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
  "psi1_0", "sigma_psi1_species",
  "psi1_herbaceous_flowers", "psi1_woody_flowers", "psi1_specialization",
  "psi1_interaction_1", "psi1_interaction_2",
  
  "gamma0", "sigma_gamma_species",
  "gamma_herbaceous_flowers", "gamma_woody_flowers", "gamma_specialization",
  "gamma_interaction_1", "gamma_interaction_2",
  
  "phi0", "sigma_phi_species",
  "phi_herbaceous_flowers", "phi_woody_flowers", "phi_specialization",
  "phi_interaction_1", "phi_interaction_2"
))


traceplot(stan_out_sim, pars = c(
  "p0", "sigma_p_species", "p_degree",
  #"sigma_p_site", 
  #"p_habitat", 
  "p_flower_abundance_any",
  "mu_p_species_date", "sigma_p_species_date", "mu_p_species_date_sq", "sigma_p_species_date_sq" 
))

### PPC's

print(stan_out_sim, digits = 3, pars = c("P_species"))

fit_summary <- rstan::summary(stan_out_sim)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
(mean_FTP <- mean(fit_summary$summary[:,1]))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)

# Evaluation of fit 
# species 1
plot(list_of_draws[,127], list_of_draws[,27], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))

abline(0, 1, lwd = 2, col = "black")

# discrepancy in the actual data too low (or discrepancy in replicate data too high)?

# site 2
plot(list_of_draws[,], list_of_draws[,], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

