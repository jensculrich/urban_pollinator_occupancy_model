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
  int<lower=0> n_years;
  int<lower=0> n_visits;
  int<lower=0,upper=1> V[n_species, n_sites, n_years, n_visits];
   
  // covariate data
  int<lower=0, upper=1> habitat_type[n_sites]; // categorical habitat type (0 or 1)
  real date_scaled[n_sites, n_years, n_visits]; // scaled day of year on which each visit was conducted

} // end data

parameters {
  
  real<lower=0, upper=1> psi1;
  real<lower=0,upper=1> gamma;
  real phi0;
  vector[n_species] phi_species;
  real sigma_phi_species;
  real phi_habitat;
  
  real p0;
  vector[n_species] p_species;
  real sigma_p_species;
  real p_habitat;
  real p_date;
  real p_date_sq;

} // end parameters

transformed parameters {
   
  // logit scale phi
  real phi[n_species, n_sites, n_years]; // odds of persistence
   
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
  
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_species[species[i]] + // a species specific persistence intercept
          phi_habitat * habitat_type[j] // a spatial persistence effect
          ); // end phi[j,k]
             
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species 
     
   
  // construct an occurrence array
  real psi[n_species, n_sites, n_years];
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        
        if(k < 2){ // define initial state
          psi[i,j,k] = psi1; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k] + (1 - psi[i,j,k-1]) * gamma; 
        } // end if/else
        
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species

  // construct an detection array
  real p[n_species, n_sites, n_years, n_visits];  // odds of detection
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
        for(l in 1:n_visits){ // loop across all visits
          
          p[i,j,k,l] = inv_logit( // probability (0-1) of detection is equal to..
            p_species[species[i]] + // a species-specific intercept
            p_habitat * habitat_type[j] + // a spatial detection effect
            // should make date species-specific
            p_date * date_scaled[j,k,l] + // a phenological detection effect (peak)
            p_date_sq * (date_scaled[j,k,l])^2 // a phenological detection effect (decay)
            ); // end p[j,k,l]
             
        } // end loop across all visits
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species 
  
   
} // end transformed parameters

model {
  
  // PRIORS
  // occupancy
  psi1 ~ uniform(0,1);
  gamma ~ uniform(0,1);
  phi0 ~ normal(0,2);
  phi_species ~ normal(phi0, sigma_phi_species); // species-specific intercepts (centered on global)
  sigma_phi_species ~ normal(0, 1); // variation in species-specific intercepts
  phi_habitat ~ normal(0,2);
  
  // detection
  p0 ~ normal(0,2); // global intercept
  p_species ~ normal(p0, sigma_p_species); // species-specific intercepts (centered on global)
  sigma_p_species ~ normal(0, 1); // variation in species-specific intercepts
  p_habitat ~ normal(0,2);
  p_date ~ normal(0,2);
  p_date_sq ~ normal(0,2);
  
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
          
      } // end loop across all intervals
    } // end loop across all sites   
  } // end loop across all species

} // end model

generated quantities{
  
  // posterior predictive check (Freeman-Tukey posterior pred check, binned by species)
  
  // estimate expected values (occurrence is a partially latent variable so when we don't observe the species, we don't actually know the expected values)
  // create replicated data using the paramter estimates and stochastic process defined by the model
  // gather real detecions
  
  // test the rate at which the number of detections in the real data versus the repped data 
  // are closer to the expected values. The FTP is the rate at which the real data are closer.
  
  int Z[n_species, n_sites, n_years];
  
  int z_rep[n_species, n_sites, n_years];
  int y_rep[n_species, n_sites, n_years, n_visits]; // repd detections

  real eval[n_species, n_sites,n_years,n_visits]; // expected values

  real T_rep[n_sites]; // Freeman-Tukey distance from eval (species bin)
  real T_obs[n_sites]; // Freeman-Tukey distance from eval (species bin)

  real P_site[n_sites]; // P-value by species

  // Initialize T_rep and T_obs and P-values
  for(j in 1:n_sites){
    
    T_rep[j] = 0;
    T_obs[j] = 0;
    
    P_site[j] = 0;

  }
      
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
      
          // if occupancy state is certain then the expected occupancy is 1
          if(sum(V[i,j,k,1:n_visits]) > 0) {
          
            Z[i,j,k] = 1;
          
          // else the site could be occupied or not
          } else {
            
            // occupancy but never observed by either dataset
            real ulo = inv_logit(psi[i,j,k]) * 
              // ((1 - inv_logit(p[j,k]))^n_visits) // for no visit heterogeneity in p
              (1 - inv_logit(p[i,j,k,1])) * (1 - inv_logit(p[i,j,k,2])) *
              (1 - inv_logit(p[i,j,k,3])) * (1 - inv_logit(p[i,j,k,4])) *
              (1 - inv_logit(p[i,j,k,5])) * (1 - inv_logit(p[i,j,k,6]));
            // non-occupancy
            real uln = (1 - inv_logit(psi[i,j,k]));
            
            // outcome of occupancy given the number of ways to get either outcome
            // higher psi values will yield a higher probability of success in the bernoulli_rng()
            // higher detection rates will yield a lower probability of success (odds are low that it was there but you never saw it)
            Z[i,j,k] = bernoulli_rng(ulo / (ulo + uln));
            
          } // end else uncertain occupancy state
        
      } // end loop across years
    } // end loop across sites
  } // end loop across species
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
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
          // Binned by species
          T_rep[j] = T_rep[j] + (sqrt(y_rep[i,j,k,l]) - 
            sqrt(eval[i,j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs[j] = T_obs[j] + (sqrt(V[i,j,k,l]) - 
            sqrt(eval[i,j,k,l]))^2;
          
        } // end loop across visits
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
  // bin by sites
  for(j in 1:n_sites) { // loop across all years
    
    // if the discrepancy is lower for the real data for the site
    // versus the replicated data
    if(T_obs[j] < T_rep[j]){
      
      // then increase site P by 1      
      P_site[j] = P_site[j] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
  }

} // end generated quantities
