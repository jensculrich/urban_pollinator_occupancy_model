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
  
  int<lower=0> n_sites;
  int<lower=0> n_years;
  int<lower=0> n_visits;
  int<lower=0,upper=1> V[n_sites, n_years, n_visits];
   
  // covariate data
  int<lower=0, upper=1> habitat_type[n_sites]; // categorical habitat type (0 or 1)
  real date_scaled[n_sites, n_years, n_visits]; // scaled day of year on which each visit was conducted

}

parameters {
  
  real p0;
  real p_habitat;
  real p_date;
  real p_date_sq;
  real<lower=0,upper=1> gamma;
  real<lower=0,upper=1> phi;
  real<lower=0, upper=1> psi1;
   
}

transformed parameters {
   
   // construct an occurrence matrix
   matrix[n_sites, n_years] psi;
   
   for (j in 1:n_sites){
     for (k in 1:n_years){
       if (k < 2){ // define initial state
          psi[j, k] = psi1; 
       } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          psi[j, k] = psi[j, k-1] * phi + (1 - psi[j, k-1]) * gamma; 
       } 
     }
   }
   
  real p[n_sites, n_years, n_visits];  // odds of detection
  
  for(j in 1:n_sites){    // loop across all sites
    for(k in 1:n_years){ // loop across all intervals
      for(l in 1:n_visits){ // loop across all visits
        
        p[j,k,l] = inv_logit( // the inverse of the log odds of detection is equal to..
          p0 + # an intercept
          p_habitat * habitat_type[j] + # a spatial detection effect
          p_date * date_scaled[j,k,l] + # a spatial detection effect
          p_date_sq * (date_scaled[j,k,l])^2  # a spatial detection effect
          ); // end p[j,k,l]
           
      } // end loop across all visits
    } // end loop across all intervals
  } // end loop across all sites
   
}

model {
  
  // priors
  psi1 ~ uniform(0,1);
  gamma ~ uniform(0,1);
  phi ~ uniform(0,1);
  
  p0 ~ normal(0,2);
  p_habitat ~ normal(0,2);
  p_date ~ normal(0,2);
  p_date_sq ~ normal(0,2);
  
  // likelihood
  for (j in 1:n_sites){
    for (k in 1:n_years){
        
        if (sum(V[j,k]) > 0){ // lp observed
          target += (log(psi[j, k]) + bernoulli_lpmf(V[j,k]|p[j,k]));
        } else { // lp unobserved (set up for 6 annual visits)
          target += (log_sum_exp(log(psi[j, k]) + log1m(p[j,k,1]) + log1m(p[j,k,2]) + 
                                                  log1m(p[j,k,3]) + log1m(p[j,k,4]) + 
                                                  log1m(p[j,k,5]) + log1m(p[j,k,6]),
                                        log1m(psi[j, k])));
        } // end if/else
        
      //}
    }
  }
}

generated quantities{
  
  // posterior predictive check (Freeman-Tukey posterior pred check, binned by species)
  
  // estimate expected values (occurrence is a partially latent variable so when we don't observe the species, we don't actually know the expected values)
  // create replicated data using the paramter estimates and stochastic process defined by the model
  // gather real detecions
  
  // test the rate at which the number of detections in the real data versus the repped data 
  // are closer to the expected values. The FTP is the rate at which the real data are closer.
  
  int Z[n_sites, n_years];
  
  int z_rep[n_sites, n_years];
  int y_rep[n_sites, n_years, n_visits]; // repd detections

  real eval[n_sites,n_years,n_visits]; // expected values

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
  //for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
      
          // if occupancy state is certain then the expected occupancy is 1
          if(sum(V[j,k,1:n_visits]) > 0) {
          
            Z[j,k] = 1;
          
          // else the site could be occupied or not
          } else {
            
            // occupancy but never observed by either dataset
            real ulo = inv_logit(psi[j,k]) * 
              // ((1 - inv_logit(p[j,k]))^n_visits) // for no visit heterogeneity in p
              (1 - inv_logit(p[j,k,1])) * (1 - inv_logit(p[j,k,2])) *
              (1 - inv_logit(p[j,k,3])) * (1 - inv_logit(p[j,k,4])) *
              (1 - inv_logit(p[j,k,5])) * (1 - inv_logit(p[j,k,6]));
            // non-occupancy
            real uln = (1 - inv_logit(psi[j,k]));
            
            // outcome of occupancy given the number of ways to get either outcome
            // higher psi values will yield a higher probability of success in the bernoulli_rng()
            // higher detection rates will yield a lower probability of success
            Z[j,k] = bernoulli_rng(ulo / (ulo + uln));
            
          } // end else uncertain occupancy state
        
      } // end loop across years
    } // end loop across sites
  //} // end loop across species
      
  // generating posterior predictive distribution
  // Predict Z at sites
  //for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
        for(l in 1:n_visits){
          
          // expected detections
          eval[j,k,l] = Z[j,k] * 
            bernoulli_logit_rng(p[j,k,l]);
          
          // occupancy in replicated data
          // should evaluate to zero if the site is not in range
          z_rep[j,k] = bernoulli_logit_rng(psi[j,k]); 

          // detections in replicated data
          y_rep[j,k,l] = z_rep[j,k] * bernoulli_logit_rng(p[j,k,l]);

          // Compute fit statistic (Tukey-Freeman) for replicate data
          // Binned by species
          T_rep[j] = T_rep[j] + (sqrt(y_rep[j,k,l]) - 
            sqrt(eval[j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs[j] = T_obs[j] + (sqrt(V[j,k,l]) - 
            sqrt(eval[j,k,l]))^2;
          
        } // end loop across visits
      } // end loop across years
    } // end loop across sites
  //} // end loop across species
  
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

}
