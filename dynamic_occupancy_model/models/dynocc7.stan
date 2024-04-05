// dynocc model 7 - dynocc 5 but with 
// categorical habitat covariate for meadows 
// (to see whether meadows after accounting for woody promote occurrence)

// dynamic multi-species occupancy model for pollinators in urban parks
// jcu, started nov 23, 2023.
// see work done here for some ideas on translating JAGS models into STAN:
// https://www.r-bloggers.com/2014/11/dynamic-occupancy-models-in-stan/

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

// let's try fitting this model with no partial-pooling of species effects
// (species all get their own intercept terms)

data {
  
  int<lower=0> n_species; // number of species
  int<lower=1> species[n_species]; // vector of species identities
  int<lower=0> n_sites; // number of sites
  int<lower=1> sites[n_sites]; // vector of site identities
  int<lower=0> n_years; // total years
  int<lower=0> n_years_minus1;
  int<lower=1> years[n_years_minus1]; // vector of year transition indexes (index 1 == transition between years 1 and 2)
  int<lower=1> years_full[n_years]; 
  int<lower=0> n_visits; // visits per year
  int<lower=0,upper=1> V[n_species, n_sites, n_years, n_visits]; // binary detection / non detection
  
  // covariate data
  int<lower=0, upper=1> habitat_type[n_sites]; // categorical habitat type (0 or 1)
  real date_scaled[n_sites, n_years, n_visits]; // scaled day of year on which each visit was conducted
  vector[n_species] d; // species specialization metric ('Bluthgen's d', but can be reset to 'degree' based on data prep)
  vector[n_species] degree; // total unique interactions made by species
  real herbaceous_flowers_scaled[n_sites, n_years]; // herbaceous flower abundance in quadrats
  real woody_flowers_scaled[n_sites, n_years]; // woody flower abundance in survey area
  real flowers_any_by_survey[n_sites, n_years, n_visits]; // flower abundance per survey visit
  
} // end data

parameters {

  // initial state
  real psi1_0; // intercept
  vector[n_species] psi1_species_raw;  
  real<lower=0> sigma_psi1_species;
  real psi1_herbaceous_flowers; // effect of herbaceous flowers
  real psi1_woody_flowers; // effect of woody flowers
  real psi1_specialization; // effect of species specialization
  real psi1_interaction_1; // interaction herbaceous flowers and specialization
  real psi1_interaction_2; // interaction woody flowers and specialization
  
  // colonization
  real gamma0;
  vector[n_species] gamma_species_raw;
  real<lower=0> sigma_gamma_species;
  real gamma_herbaceous_flowers;
  real gamma_woody_flowers;
  real gamma_specialization;
  real gamma_interaction_1;
  real gamma_interaction_2;
  vector[n_years_minus1] gamma_year;
  
  // persistence
  real phi0;
  vector[n_species] phi_species_raw;
  real<lower=0> sigma_phi_species;
  real phi_herbaceous_flowers;
  real phi_woody_flowers;
  real phi_specialization;
  real phi_interaction_1;
  real phi_interaction_2;
  vector[n_years_minus1] phi_year;

  // detection
  real p0; // intercept
  vector[n_species] p_species_raw;
  real<lower=0> sigma_p_species;
  real p_specialization; // effect of species specialization
  vector[n_species] p_date; // phenology peak
  real mu_p_species_date; // community mean
  real<lower=0> sigma_p_species_date; // variation
  vector[n_species] p_date_sq; // decay pattern of phenology
  real mu_p_species_date_sq; // variation
  real<lower=0> sigma_p_species_date_sq; // community mean
  real p_flower_abundance_any; // effect of survey specific flower abundance
  vector[n_years] p_year;

} // end parameters

transformed parameters {

  // logit scale psi1, gamma, phi
  real psi1[n_species, n_sites]; // odds of occurrence year 1
  real gamma[n_species, n_sites, n_years_minus1]; // odds of colonization
  real phi[n_species, n_sites, n_years_minus1]; // odds of persistence
  real p[n_species, n_sites, n_years, n_visits];  // odds of detection
  
  vector[n_species] psi1_species;
  vector[n_species] gamma_species;
  vector[n_species] phi_species;
  vector[n_species] p_species;
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_species = psi1_0 + sigma_psi1_species * psi1_species_raw;
  gamma_species = gamma0 + sigma_gamma_species * gamma_species_raw;
  phi_species = phi0 + sigma_phi_species * phi_species_raw;
  p_species = p0 + sigma_p_species * p_species_raw;
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years_minus1){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          psi1_species[species[i]] + // a species specific intercept
          psi1_specialization * d[i] +
          psi1_herbaceous_flowers * habitat_type[j] + 
          (psi1_interaction_1 * d[i] * habitat_type[j]) +
          psi1_woody_flowers * woody_flowers_scaled[j,1] + 
          (psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
          ); // end phi[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma_species[species[i]] + // a species specific intercept
          gamma_year[years[k]] + // effect of first year transition // index of years 1 == time between year 1 and 2
          gamma_specialization * d[i] +
          // colonization depending on incoming year's site conditions (k+1)
          gamma_herbaceous_flowers * habitat_type[j] + 
          (gamma_interaction_1 * d[i] * habitat_type[j]) +
          gamma_woody_flowers * woody_flowers_scaled[j,k+1] + 
          (gamma_interaction_2 * d[i] * woody_flowers_scaled[j,k+1])
          ); // end phi[j,k]
                
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_species[species[i]] + // a species specific intercept
          phi_year[years[k]] + // the effect of the first year transition // index 1 == year 2
          phi_specialization * d[i] +
          // persistence depending on outgoing year's site conditions (k)
          phi_herbaceous_flowers * habitat_type[j] + 
          (phi_interaction_1 * d[i] * habitat_type[j]) +
          phi_woody_flowers * woody_flowers_scaled[j,k] + 
          (phi_interaction_2 * d[i] * woody_flowers_scaled[j,k]) 
          ); // end phi[j,k]
          
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
  
  // have to do another loop because phi and gamma have a shorter k index length than p
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_visits){ // loop across all visits

          p[i,j,k,l] = inv_logit( // probability (0-1) of detection is equal to..
            p_species[species[i]] + // a species specific intercept
            p_year[years_full[k]] + 
            p_specialization * degree[i] +
            p_date[species[i]] * date_scaled[j,k,l] + // a species-specific phenological detection effect (peak)
            p_date_sq[species[i]] * (date_scaled[j,k,l])^2 + // a species-specific phenological detection effect (decay)
            p_flower_abundance_any * flowers_any_by_survey[j,k,l]
            ); // end p[j,k,l]
            
        } // end loop across all visits
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
   
  // construct an occurrence array
  real psi[n_species, n_sites, n_years];
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        
        if(k < 2){ // define initial state
          psi[i,j,k] = psi1[i,j]; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          // reduce 1 from k for phi and gamma because there are n_years - 1 transitions 
          // and so there are only n_years - 1 speciesXsite "stacks" of phi and gamma
          // but phi[,,k-1] for k = 2 will actually consider the effects of 
          // e.g. flower abundacnce in year 2 (since year 2 is the first year we estimate phi)
          psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k-1] + (1 - psi[i,j,k-1]) * gamma[i,j,k-1]; 
        } // end if/else
        
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 2); // persistence intercept
  psi1_species_raw ~ std_normal();
  sigma_psi1_species ~ normal(0, 2);
  psi1_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  psi1_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  psi1_specialization ~ normal(0, 2);
  psi1_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  psi1_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  
  // colonization
  gamma0 ~ normal(0, 1); // persistence intercept
  gamma_species_raw ~ std_normal();
  sigma_gamma_species ~ normal(0, 1); 
  gamma_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  gamma_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  gamma_specialization ~ normal(0, 2);
  gamma_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  gamma_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  gamma_year ~ normal(0, 0.25); // year effects
  
  // persistence
  phi0 ~ normal(0, 1); // global intercept
  phi_species_raw ~ std_normal();
  sigma_phi_species ~ normal(0, 1);
  phi_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  phi_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  phi_specialization ~ normal(0, 2);
  phi_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  phi_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  phi_year ~ normal(0, 0.25);
  
  // detection
  p0 ~ normal(0, 2); // global intercept
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 2);
  p_specialization ~ normal(0, 2); // effect of specialization on intercept
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0, 2); // mean
  sigma_p_species_date ~ cauchy(0, 2); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0, 2); // mean
  sigma_p_species_date_sq ~ cauchy(0, 2); // variation
  p_flower_abundance_any ~ normal(0, 2); // effect of survey visit flower abundance on detection
  p_year ~ normal(0, 0.25);
  
  // LIKELIHOOD
  for(i in 1:n_species){
    for (j in 1:n_sites){
      for (k in 1:n_years){
          
          if (sum(V[i,j,k]) > 0){ // lp observed 
            // detection on each visit given detection rate on each visit
            target += (log(psi[i,j,k]) + 
                       + bernoulli_lpmf(V[i,j,k]|p[i,j,k]));
          } else { // lp unobserved (set up for 6 annual visits)
            // marginal likelihood of...
            // non-detection on each visit given detection rate on each visit given occurrence
            target += (log_sum_exp(log(psi[i,j,k]) + log1m(p[i,j,k,1]) + log1m(p[i,j,k,2]) + 
                                                    log1m(p[i,j,k,3]) + log1m(p[i,j,k,4]) + 
                                                    log1m(p[i,j,k,5]) + log1m(p[i,j,k,6]),
            // or just simple no occurrence
                                          log1m(psi[i,j,k])));
          } // end if/else
  
      } // end loop across all years
    } // end loop across all sites   
  } // end loop across all species

} // end model

generated quantities{
  
  // Could also add: 
  // effects on init occupancy, colonization, and persistence across a range of specialization bins
  
  // Diversity estimation
  // number of species at each site in each year
  int z_simmed[n_species, n_sites, n_years]; // simulate occurrence
  int species_richness[n_sites, n_years]; // site/year species richness
  real avg_species_richness_enhanced[n_years]; // average across sites
  real avg_species_richness_control[n_years]; // average across sites
  real increase_richness_enhanced[n_years]; // difference in average species richness in each year
  
  for(i in 1:n_species){
   for(j in 1:n_sites){
     for(k in 1:n_years){
          z_simmed[i,j,k] = bernoulli_rng(psi[i,j,k]); 
      }    
    }
  }
  
  // Initialize avg_species_richness
  for(k in 1:n_years){
    avg_species_richness_control[k] = 0;
    avg_species_richness_enhanced[k] = 0;
  }
  
  // summed number of species occurring in site/year
  for(j in 1:n_sites){
    for(k in 1:n_years){
      
      // first calc site/year specific species richness
      species_richness[j,k] = sum(z_simmed[,j,k]);
      
      // use habitat type as a switch for whether or not you will add it to the species richness group
      avg_species_richness_control[k] = avg_species_richness_control[k] + ((1 - habitat_type[j]) *  species_richness[j,k]);
      avg_species_richness_enhanced[k] = avg_species_richness_enhanced[k] + (habitat_type[j] *  species_richness[j,k]);
      
    }
  }
  
  for(k in 1:n_years){
    // average out species richness by number of sites (1/2 of total sites were in each category)
    avg_species_richness_control[k] = avg_species_richness_control[k] / (n_sites / 2.0);
    avg_species_richness_enhanced[k] = avg_species_richness_enhanced[k] / (n_sites / 2.0);
  }
  
  for(k in 1:n_years){
    increase_richness_enhanced[k] = avg_species_richness_enhanced[k] - avg_species_richness_control[k];
  }
  
  //
  // posterior predictive check (number of detections, binned by species)
  //
  //int<lower=0,upper=1> V_rep[n_species, n_sites, n_years, n_visits]; // repd detections given repd occurrence
  int<lower=0> W_species_rep[n_species]; // sum of simulated detections
  
  // initialize at 0
  for(i in 1:n_species){
    W_species_rep[i] = 0;
  }
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_visits){
          
          // detections in replicated data (us z_simmed from above)
          W_species_rep[i] = W_species_rep[i] + 
            (z_simmed[i,j,k] * bernoulli_rng(p[i,j,k,l]));
           
        } // end loop across visits
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
} // end generated quantities
