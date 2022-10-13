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
### Variable values for data simulation
## study dimensions
n_species = 25 ## number of species
n_sites = 25 ## number of sites
n_intervals = 4 ## number of occupancy intervals (years)
n_visits = 7 ## number of samples per year

## occupancy
mu_psi_0 = -.5 # global occupancy intercept

# sigma_psi_species = 0.5 # species specific occupancy intercept

# mu_psi_species_habitat # community mean effect of habitat on occupancy
# sigma_psi_species_habitat # variation across species in the effect of habitat on occupancy

# psi_interval # fixed effect of year on occupancy (to check whether it varies between years)

## detection
mu_p_0 = -1 # global detection intercept

# mu_p_species_date # community mean effect of day of year on detection
# sigma_p_species_date # variation across species in the effect of day of year on detection

# mu_p_species_date_sq # community mean effect of (day of year)^2 on detection (enabling a seasonal peak)
# sigma_p_species_date_sq # variation across species in the effect of day of year on detection

# mu_p_species_habitat = 0.25 # community mean effect of habitat on detection
# sigma_p_species_habitat = 0.3 # variation across species in the effect of habitat on detection

## expit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## Create covariate data

## habitat
# should be 0 (mowed) or 1 (meadow)
# should have one value per site

## day of year
# should have a value for each site*year*visit

# Interval values: numeric vector (will act as covariate data for psi.interval)
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)


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
        mu_psi_0 # + # a baseline intercept
          # psi_species[species] + # a species specific intercept
          # mu_psi_species_habitat[species]*habitat[site] # a species specific effect of habitat
          # psi_interval[species]*intervals[interval] # a species specific temporal change
      )
      
      for(visit in 1:n_visits) { # for each visit
        
        p_matrix[species, site, interval, visit] <- ilogit( # detection is equal to 
          mu_p_0 # + # a baseline intercept
            # p_species[species] + # a species specific intercept
            # mu_p_species_habitat[species]*habitat[site]
            # mu_p_species_date[species]*date[visit] + # a spatiotemporally specific intercept
            # mu_p_species_date_sq[species]*(date[visit])^2 + # a spatiotemporally specific intercept
        )
        
      } # for each visit
    } # for each species
  } # for each interval
} # for each site

# preview the psi and p arrays
head(psi_matrix[1:n_species, 1:n_sites,1])
head(p_matrix[1:n_species, 1:n_sites,1,1])

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

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- V # detection data
n_species <- n_species # number of species
n_sites <- n_sites # number of sites
n_intervals <- n_intervals # number of surveys 
n_visits <- n_visits

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

stan_data <- c("V", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites") 

# Parameters monitored
params <- c("mu_psi_0",
            #"psi_species",
            # "sigma_psi_species",
            # "psi_interval",
            # "mu_psi_interval",
            # "sigma_psi_interval",
            "mu_p_0"# ,
            # "p_species",
            # "sigma_p_species",
            # "p_site",
            # "sigma_p_site",
            # "p_interval"
)


# MCMC settings
n_iterations <- 500
n_thin <- 1
n_burnin <- 250
n_chains <- 4
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       #sigma_psi_species = runif(1, 0, 1),
       #mu_psi_interval = runif(1, -1, 1),
       #sigma_psi_interval = runif(1, 0, 1),
       mu_p_0 = runif(1, -1, 1)#,
       #sigma_p_species = runif(1, 0, 1),
       #sigma_p_site = runif(1, 0, 1),
       #p_interval = runif(1, -1, 1)
       
  )
)

#targets <- as.data.frame(cbind(param_names, parameters))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./simulation/model0.stan"

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

