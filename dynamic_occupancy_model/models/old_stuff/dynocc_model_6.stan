// dynamic multi-species occupancy model for pollinators in urban parks
// jcu, started nov 23, 2023.
// see work done here for some ideas on translating JAGS models into STAN:
// https://www.r-bloggers.com/2014/11/dynamic-occupancy-models-in-stan/

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

data {
  
  int<lower=0> n_species;
  int<lower=1> species[n_species]; // vector of species identities
  int<lower=0> n_sites;
  int<lower=1> sites[n_sites]; // vector of site identities
  int<lower=0> n_years;
  int<lower=0> n_years_minus1;
  int<lower=0> years[n_years_minus1]; // vector of year identities (no final year because we don't return after last transition)
  int<lower=0> n_visits;
  int<lower=0,upper=1> V[n_species, n_sites, n_years, n_visits];
  
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
  
  // correlated species random effects
  vector<lower=0>[4] sigma_species; // SDs for random effects for persistence
  cholesky_factor_corr[4] L_species; // Correlation matrix for random intercepts and slopes for persistence
  matrix[4,n_species] z_species; // Random effects for persistence//real delta1_phi0;
  
  // initial state
  real psi1_0; 
  real psi1_herbaceous_flowers;
  real psi1_woody_flowers;
  real psi1_specialization;
  real psi1_interaction_1;
  real psi1_interaction_2;
  
  // colonization
  real gamma0; 
  real gamma_herbaceous_flowers;
  real gamma_woody_flowers;
  real gamma_specialization;
  real gamma_interaction_1;
  real gamma_interaction_2;
  vector[n_years_minus1] gamma_year;
  
  // persistence
  real phi0; 
  real phi_herbaceous_flowers;
  real phi_woody_flowers;
  real phi_specialization;
  real phi_interaction_1;
  real phi_interaction_2;
  vector[n_years_minus1] phi_year;

  // detection
  real p0;
  real p_specialization;
  vector[n_species] p_date;
  real mu_p_species_date;
  real<lower=0> sigma_p_species_date;
  vector[n_species] p_date_sq;
  real<lower=0> sigma_p_species_date_sq;
  real mu_p_species_date_sq;
  real p_flower_abundance_any;

} // end parameters

transformed parameters {

  // logit scale psi1, gamma, phi
  real psi1[n_species, n_sites]; // odds of occurrence year 1
  real gamma[n_species, n_sites, n_years]; // odds of colonization
  real phi[n_species, n_sites, n_years]; // odds of persistence
  
  matrix[n_species,4] species_effects;
  species_effects = (diag_pre_multiply(sigma_species,L_species)*z_species)'; // the dash indicates transposition
  
   
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          psi1_0 +
          species_effects[species[i],1] + // a species specific intercept
          psi1_herbaceous_flowers * herbaceous_flowers_scaled[j,1] + 
          psi1_specialization * d[i] +
          (psi1_interaction_1 * d[i] * herbaceous_flowers_scaled[j,1]) +
          psi1_woody_flowers * woody_flowers_scaled[j,1] + 
          (psi1_interaction_2 * d[i] * woody_flowers_scaled[j,1])
          ); // end phi[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma0 +
          species_effects[species[i],2] + // a species specific intercept
          gamma_herbaceous_flowers * herbaceous_flowers_scaled[j,k] + 
          gamma_specialization * d[i] +
          (gamma_interaction_1 * d[i] * herbaceous_flowers_scaled[j,k]) +
          gamma_woody_flowers * woody_flowers_scaled[j,k] + 
          (gamma_interaction_2 * d[i] * woody_flowers_scaled[j,k]) +
          gamma_year[years[k]]
          ); // end phi[j,k]
        
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          phi0 +
          species_effects[species[i],3] + // a species specific intercept
          phi_herbaceous_flowers * herbaceous_flowers_scaled[j,k] + 
          phi_specialization * d[i] +
          (phi_interaction_1 * d[i] * herbaceous_flowers_scaled[j,k]) +
          phi_woody_flowers * woody_flowers_scaled[j,k] + 
          (phi_interaction_2 * d[i] * woody_flowers_scaled[j,k]) +
          phi_year[years[k]]
          ); // end phi[j,k]
             
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
          psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k] + (1 - psi[i,j,k-1]) * gamma[i,j,k]; 
        } // end if/else
        
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species

  // construct an detection array
  real p[n_species, n_sites, n_years, n_visits];  // odds of detection
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_visits){ // loop across all visits
          
          p[i,j,k,l] = inv_logit( // probability (0-1) of detection is equal to..
            p0 +
            species_effects[species[i],4] + + // a species specific intercept
            p_specialization * degree[i] +
            p_date[species[i]] * date_scaled[j,k,l] + // a species-specific phenological detection effect (peak)
            p_date_sq[species[i]] * (date_scaled[j,k,l])^2 + // a species-specific phenological detection effect (decay)
            p_flower_abundance_any * flowers_any_by_survey[j,k,l]
            ); // end p[j,k,l]
             
        } // end loop across all visits
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
  
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // correlated species random effects
  L_species ~ lkj_corr_cholesky(2); // Ranef prior
  z_species[1,] ~ normal(0, 1);
  z_species[2,] ~ normal(0, 1);
  z_species[3,] ~ normal(0, 1);
  z_species[4,] ~ normal(0, 1);
  sigma_species[1] ~ normal(0, 1);
  sigma_species[2] ~ normal(0, 1);
  sigma_species[3] ~ normal(0, 1);
  sigma_species[4] ~ normal(0, 1);
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 2); // persistence intercept
  psi1_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  psi1_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  psi1_specialization ~ normal(0, 2);
  psi1_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  psi1_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  
  // colonization
  gamma0 ~ normal(0, 1); // persistence intercept
  gamma_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  gamma_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  gamma_specialization ~ normal(0, 2);
  gamma_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  gamma_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  gamma_year ~ normal(0, 1); // year effects
  
  // persistence
  phi0 ~ normal(0, 1); // persistence intercept
  phi_herbaceous_flowers ~ normal(0, 2); // effect of habitat on colonization
  phi_woody_flowers ~ normal(0, 2); // effect of habitat on colonization
  phi_specialization ~ normal(0, 2);
  phi_interaction_1 ~ normal(0, 1); // baseline effect of habitat 
  phi_interaction_2 ~ normal(0, 1); // effect of specialization on response to habitat
  phi_year ~ normal(0, 1);
  
  // detection
  p0 ~ normal(0, 2); // global intercept
  p_specialization ~ normal(0, 2); // effect of specialization on intercept
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0, 2); // mean
  sigma_p_species_date ~ normal(0, 1); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0, 2); // mean
  sigma_p_species_date_sq ~ normal(0, 1); // variation
  p_flower_abundance_any ~ normal(0, 2); // effect of survey visit flower abundance on detection
  
  // LIKELIHOOD
  for(i in 1:n_species){
    for (j in 1:n_sites){
      for (k in 1:n_years){
          
          if (sum(V[i,j,k]) > 0){ // lp observed
            target += (log(psi[i,j,k]) + bernoulli_lpmf(V[i,j,k]|p[i,j,k]));
          } else { // lp unobserved (set up for 6 annual visits)
            target += (log_sum_exp(log(psi[i,j,k]) + log1m(p[i,j,k,1]) + log1m(p[i,j,k,2]) + 
                                                    log1m(p[i,j,k,3]) + log1m(p[i,j,k,4]) + 
                                                    log1m(p[i,j,k,5]) + log1m(p[i,j,k,6]),
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
  
  // equilibrium occupancy rate (for average species)  
  //real psi_eq_habitat0; // control sites
  //real psi_eq_habitat1; // enhanced sites
  
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
  
  //psi_eq_habitat0 = inv_logit(gamma0) / (inv_logit(gamma0)+(1-inv_logit(delta0_phi0))); // equilibrium occupancy rate 
  //psi_eq_habitat1 = inv_logit(gamma0 + mu_gamma_habitat) / (inv_logit(gamma0 + mu_gamma_habitat)+(1-inv_logit(delta0_phi0 + delta0_phi_habitat))); // equilibrium occupancy rate
  
  // annual turnover rates
  // need to figure out how to fix this for spatially and species varying psi and gamma
  //real turnover_control[n_years-1]; // average across sites
  //real turnover_enhanced[n_years-1]; // average across sites
  
  //for(k in 2:n_years){
  //  turnover_control[k-1] = (1 - psi[,,k-1]) *  gamma[,,k-1])/psi[,,k];
  //  turnover_enhanced[k-1] = (1 - psi[,,k-1]) * gamma[,,k-1])/psi[,,k];
  //}
  
  //
  // posterior predictive check (Freeman-Tukey posterior pred check, binned by species)
  //
  // estimate expected values (occurrence is a partially latent variable so when we don't observe the species, we don't actually know the expected values)
  // create replicated data using the paramter estimates and stochastic process defined by the model
  // gather real detecions
  
  // test the rate at which the number of detections in the real data versus the repped data 
  // are closer to the expected values. The FTP is the rate at which the real data are closer.
  
  int Z[n_species, n_sites, n_years];
  
  int z_rep[n_species, n_sites, n_years]; // repd occurrence
  int y_rep[n_species, n_sites, n_years, n_visits]; // repd detections given repd occurrence

  real eval[n_species,n_sites,n_years,n_visits]; // expected values

  real T_rep[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs[n_species]; // Freeman-Tukey distance from eval (species bin)

  real P_species[n_species]; // P-value by species

  // Initialize T_rep and T_obs and P-values
  for(i in 1:n_species){
    
    T_rep[i] = 0;
    T_obs[i] = 0;
    
    P_species[i] = 0;

  }
      
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all years
      
          // if occupancy state is certain then the expected occupancy is 1
          if(sum(V[i, j, k, 1:n_visits]) > 0) {
          
            Z[i,j,k] = 1;
          
          // else the site could be occupied or not
          } else {
            
            //Z[i,j,k] = bernoulli_logit_rng(psi[i,j,k]);
            
            // occupancy but never observed by either dataset
            real ulo = inv_logit(psi[i,j,k]) * 
              //((1 - inv_logit(p[j,k]))^n_visits) // for no visit heterogeneity in p
              (1 - inv_logit(p[i,j,k,1])) * (1 - inv_logit(p[i,j,k,2])) *
              (1 - inv_logit(p[i,j,k,3])) * (1 - inv_logit(p[i,j,k,4])) *
              (1 - inv_logit(p[i,j,k,5])) * (1 - inv_logit(p[i,j,k,6]));
            // non-occupancy
            real uln = (1 - inv_logit(psi[i,j,k]));
            
            // outcome of occupancy given the likelihood associated with both possibilities
            Z[i,j,k] = bernoulli_rng(ulo / (ulo + uln));
            
          } // end else uncertain occupancy state
        
      } // end loop across years
    } // end loop across sites
  } // end loop across species
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_visits){
          
          // expected detections
          eval[i,j,k,l] = Z[i,j,k] * 
            bernoulli_logit_rng(p[i,j,k,l]);
          
          // occupancy in replicated data
          // should evaluate to zero if the site is not in range
          z_rep[i,j,k] = bernoulli_logit_rng(psi[i,j,k]); 
          
          // detections in replicated data
          y_rep[i,j,k,l] = z_rep[i,j,k] * bernoulli_logit_rng(p[i,j,k,l]);

          // Compute fit statistic (Tukey-Freeman) for replicate data
          // community science records
          // Binned by species
          T_rep[i] = T_rep[i] + (sqrt(y_rep[i,j,k,l]) - 
            sqrt(eval[i,j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs[i] = T_obs[i] + (sqrt(V[i,j,k,l]) - 
            sqrt(eval[i,j,k,l]))^2;
          
        } // end loop across visits
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
  // bin by species
  for(i in 1:n_species) { // loop across all species
    
    // if the discrepancy is lower for the real data for the species
    // versus the replicated data
    if(T_obs[i] < T_rep[i]){
      
      // then increase species P by 1      
      P_species[i] = P_species[i] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
  }

} // end generated quantities
