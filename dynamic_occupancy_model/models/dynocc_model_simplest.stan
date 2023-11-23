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
  
  //real<lower=0,upper=1> p;
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
       if (k < 2){ // initial state
          psi[j, k] = psi1; 
       } else { // system dynamics
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
          p_date_sq * (date_scaled[j,k,l]^2)  # a spatial detection effect
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
  //p ~ uniform(0,1);
  p0 ~ normal(0,2);
  p_habitat ~ normal(0,2);
  p_date ~ normal(0,2);
  p_date_sq ~ normal(0,2);
  
  // likelihood
  for (j in 1:n_sites){
    for (k in 1:n_years){
      for (l in 1:n_visits){
        
        if (sum(V[j, k]) > 0){
          target += (log(psi[j, k]) + bernoulli_lpmf(V[j,k,l]|p[j,k,l]));
        } else {
          target += (log_sum_exp(log(psi[j, k]) + bernoulli_lpmf(V[j,k,l]|p[j,k,l]),
                                        log(1-psi[j, k])));
        }
        
      }
    }
  }
}
