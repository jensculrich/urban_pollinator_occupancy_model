### Data simulation for pollinator populations in urban parks
# jcu; started oct 12, 2022

# multiple species will be simulated hierarchically across spatially replicated sites
# half of the sites in habitat category 0 and the other half in category 1.
# the data will be simulated to be from multiple single years (occupancy intervals)
# that are each assumed closed (no extinction/colonization in intervals), 
# and not yet linked in a dynamic form. 
# each site will be experience several temporally replicated visits within each year.
# species will be simulated such that they do not occur uniformly across the year,
# rather have a phenological window with a seasonal peak.

## --------------------------------------------------
## Operation Functions

## expit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## --------------------------------------------------
### Make Data Function
make_data <- function(## study dimensions
  n_species = 35, ## number of species
  n_sites = 30, ## number of sites
  n_intervals = 3, ## number of occupancy intervals (years)
  n_visits = 6, ## number of samples per year
  
  # introduce NAs (missed visits)?
  introduce_NAs = FALSE,
  
  ## occupancy
  mu_psi_0 = -.5, # global occupancy intercept
  
  sigma_psi_species = 0.5, # species specific occupancy intercept
  
  mu_psi_species_habitat = 1, # community mean effect of habitat on occupancy
  sigma_psi_species_habitat = 0.5, # variation across species in the effect of habitat on occupancy
  
  # psi_interval # fixed effect of year on occupancy (to check whether it varies between years)
  
  ## detection
  mu_p_0 = 0, # global detection intercept
  
  sigma_p_species = 0.5, # species specific detection intercept
  
  mu_p_species_date = 0, # community mean effect of day of year on detection
  sigma_p_species_date = 0.5, # variation across species in the effect of day of year on detection
  
  mu_p_species_date_sq = -0.5, # community mean effect of (day of year)^2 on detection (enabling a seasonal peak)
  sigma_p_species_date_sq = 0.5, # variation across species in the effect of day of year on detection
  
  mu_p_species_habitat = 0, # community mean effect of habitat on detection
  sigma_p_species_habitat = 0.5, # variation across species in the effect of habitat on detection
  
  # survey date inputs
  mean_survey_date = 180,
  sigma_survey_dates = 40) 

{

  ## --------------------------------------------------
  ## Create covariate data
  
  ## day of year
  # should have a value for each site*year*visit
  date <- array(NA, dim =c(n_sites, n_intervals, n_visits))
  
  for(site in 1:n_sites) { # for each site
    for(interval in 1:n_intervals) { # for each interval
      # create a vector of visit dates centered on the middle of the early summer
      date[site, interval, 1:n_visits] <- sort(as.integer(rnorm(
        n_visits, mean = mean_survey_date, sd = sigma_survey_dates))) 
    }
  }
  
  ## let's scale the calendar day of year by the mean date (z-score scaled)
  date_scaled <- center_scale(date) 
  
  #head(date[1,1,1:n_visits])
  #head(date_scaled[1,1,1:n_visits])
  #head(date[1,2,1:n_visits])
  #head(date_scaled[12,1,1:n_visits])
  
  ## habitat
  # should be 0 (mowed) or 1 (meadow)
  # should have one value per site
  
  habitat_category <- sort(rep(c(0, 1), times = (n_sites/2))) # half sites in each category
  
  # Interval values: numeric vector (will act as covariate data for psi.interval)
  intervals_raw <- seq(1, n_intervals, by=1)
  intervals <- intervals_raw - 1
  sites <- seq(1, n_sites, by=1)
  species <- seq(1, n_species, by=1)
  
  ## --------------------------------------------------
  ### specify species-specific occupancy probabilities
  
  ## species specific occupancy intercept (some species occur more frequently than others)
  psi_species <- rnorm(n=n_species, mean=0, sd=sigma_psi_species)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_psi_species
  
  ## effect of habitat category (0 or 1) on species occurrence (species-specific random slope)
  psi_species_habitat <- rnorm(n=n_species, mean=mu_psi_species_habitat, sd=sigma_psi_species_habitat)
  # change in each species occurrence probability from site type 0 to 1 is drawn from 
  # a distribution defined by a community mean (mu_psi_species_habitat) with 
  # species specific variation defined by sigma_psi_species_habitat
  
  ## --------------------------------------------------
  ### specify species-specific detection probabilities
  
  ## species specific occupancy intercept (some species easier to detect than others)
  p_species <- rnorm(n=n_species, mean=0, sd=sigma_p_species)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_psi_species
  
  ## effect of habitat category (0 or 1) on species detection (species-specific random slope)
  p_species_habitat <- rnorm(n=n_species, mean=mu_p_species_habitat, sd=sigma_p_species_habitat)
  # change in each species occurrence probability from site type 0 to 1 is drawn from 
  # a distribution defined by a community mean (mu_psi_species_habitat) with 
  # species specific variation defined by sigma_psi_species_habitat
  
  ## effect of date on detection (species-specific random slope)
  p_species_date <- rnorm(n=n_species, mean=mu_p_species_date, sd=sigma_p_species_date)
  # change in each species detection probability through the season is drawn from 
  # a distribution defined by a community mean (mu_p_species_date) with 
  # species specific variation defined by sigma_p_species_date
  
  ## effect of date on detection (species-specific random slope)
  p_species_date_sq <- rnorm(n=n_species, mean=mu_p_species_date_sq, sd=sigma_p_species_date_sq)
  # change in each species detection probability through the season is drawn from 
  # a distribution defined by a community mean (mu_p_species_date) with 
  # species specific variation defined by sigma_p_species_date
  
  ## --------------------------------------------------
  ## Create arrays for psi and p
  psi_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  p_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        psi_matrix[species, site, interval] <- ilogit( # occupancy is equal to
          mu_psi_0 + # a baseline intercept
            psi_species[species] + # a species specific intercept
            psi_species_habitat[species]*habitat_category[site] # a species specific effect of habitat
          # should eventually change above predictor to include a random effect of species,site (not just species) 
          # to account for multiple years worth of data (pseudoreplicates) from each site.
          # may not be necessary if the model is developed into a dynamic form linking closed seasons across years.
          # psi_interval[species]*intervals[interval] # a species specific temporal change
        )
        
        for(visit in 1:n_visits) { # for each visit
          
          p_matrix[species, site, interval, visit] <- ilogit( # detection is equal to 
            mu_p_0 + # a baseline intercept
              p_species[species] + # a species specific intercept
              p_species_habitat[species]*habitat_category[site] +
              p_species_date[species]*date_scaled[site, interval, visit] + # a spatiotemporally specific intercept
              p_species_date_sq[species]*(date_scaled[site, interval, visit])^2 # a spatiotemporally specific intercept
          )
          
        } # for each visit
      } # for each species
    } # for each interval
  } # for each site
  
  # preview the psi and p arrays
  # head(psi_matrix[1:n_species, 1:n_sites,1]) # should be higher (for the average species) 
  # in the sites coded as habitat category 1 (if input value of mu_psi_species_habitat > 0)
  
  # head(p_matrix[which.max(p_species_date),1,1,1:n_visits]) # should have high detection prob at late visits
  # head(p_matrix[which.min(p_species_date),1,1,1:n_visits]) # should have low detection prob at earlier visits
  
  
  ## --------------------------------------------------
  ## Generate presence and absence
  
  Z <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        
        # eventually here we will need to specify which(site) to restrict to actual range
        Z[species,site,interval] <- rbinom(n = 1, size = 1, 
                                           prob = psi_matrix[species,site,interval])
        
      }
    }
  }
  
  # head(Z[1,1,1])
  
  ## --------------------------------------------------
  ## Generate detection non detection data 
  
  V <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals,
                       n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          # eventually here we will need to specify which(site) to restrict to actual range
          V[species,site,interval,visit] <- Z[species,site,interval] * # occupancy state * detection prob
            rbinom(n = 1, size = 1, prob = p_matrix[species,site,interval,visit])
          
        }
      }
    }
  }
  
  if(introduce_NAs == TRUE){
    
    interval_missed = sample.int(n_intervals, 1)
    visit_missed = sample.int(n_visits, 1)
    
    V[1:n_species, 1:n_sites, interval_missed, visit_missed] <- NA
    
  }
  
  # head(V_NA[1:n_species, 1:n_sites, interval_missed, visit_missed])
  
  # Return stuff
  return(list(
    n_species = n_species, ## number of species
    n_sites = n_sites, ## number of sites
    n_intervals = n_intervals, ## number of occupancy intervals (years)
    n_visits = n_visits, ## number of samples per year
    V = V,
    date_scaled = date_scaled,
    habitat_category = habitat_category
   ))
  
}


## --------------------------------------------------
### Prepare data for model
my_data <- make_data()
str(my_data)

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 
n_visits <- my_data$n_visits

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

date_scaled <- my_data$date_scaled
habitat_category <- my_data$habitat_category

stan_data <- c("V", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites", 
               "habitat_category", "date_scaled") 

# Parameters monitored
params <- c("mu_psi_0",
            "sigma_psi_species",
            "mu_psi_species_habitat",
            "sigma_psi_species_habitat",
            "mu_p_0",
            "sigma_p_species",
            "mu_p_species_habitat",
            "sigma_p_species_habitat",
            "mu_p_species_date",
            "sigma_p_species_date",
            "mu_p_species_date_sq",
            "sigma_p_species_date_sq"
)


# MCMC settings
n_iterations <- 600
n_thin <- 1
n_burnin <- 300
n_chains <- 3
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       mu_psi_species_habitat = runif(1, -1, 1),
       sigma_psi_species_habitat = runif(1, 0, 1),
       mu_p_0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       mu_p_species_habitat = runif(1, -1, 1),
       sigma_p_species_habitat = runif(1, 0, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 1),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

#targets <- as.data.frame(cbind(param_names, parameters))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./models/model0.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

print(stan_out_sim, digits = 3)
saveRDS(stan_out_sim, "./simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim.rds")

traceplot(stan_out_sim, pars = c(
          "mu_psi_0",
          "sigma_psi_species",
          "mu_p_0", 
          "sigma_p_species",
          "mu_p_species_date",
          "sigma_p_species_date",
          "mu_p_species_date_sq",
          "sigma_p_species_date_sq"))
