// dynamic multi-species occupancy model for pollinators in urban parks
// jcu, started nov 23, 2023.
// see work done here for some ideas on translating JAGS models into STAN:
// https://github.com/stan-dev/example-models/blob/master/BPA/Ch.13/Dynocc.stan#L74
// there are some tricky things M. Joseph has done to build the model template
// will need to walk through to understand fully

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

data {
  
  int<lower=1> n_sites;  // sites within region
  //int<lower=1> sites[n_sites];  // vector of sites
  
  int<lower=1> n_years; // number of years surveyed
  int<lower=1> n_visits; // number of surveys per year
  
  int<lower=0> V[n_sites, n_years, n_visits]; // detection / non-detection data
  
  int<lower=0, upper=1> habitat_type[n_sites]; // categorical habitat type (0 or 1)
  
} // end data

transformed data {
  
  array[n_sites, n_years] int<lower=1, upper=n_visits + 1> sum_y; // sum + 1
  
  // calculate total years with transitions
  int ny_minus_1 = n_years - 1;
  
  // sum_y is the number of detections (+1) for a site X year
  for (j in 1 : n_sites) {
    for (k in 1 : n_years) {
      sum_y[j, k] = sum(V[j, k, 1 : n_visits]) + 1;
    }
  }
  
} // end  transformed data

parameters {
  
  real<lower=0, upper=1> psi1; // Occupancy probability at k=1
  
  vector<lower=0, upper=1>[nyear] p; // Detection probability
  real p0; // Detection probability (intercept for each year)
  real p_habitat_type; // Detection probability (effect of habitat type)
  
  array[2, n_years - 1] simplex[2] ps; // Transition probability
  // This is equivalent to the following:
  //  ps[1, k, 1] = phi[k];
  //  ps[1, k, 2] = 1.0 - phi[k];
  //  ps[2, k, 1] = gamma[k];
  //  ps[2, k, 2] = 1.0 - gamma[k];

}

transformed parameters {
  
  real p[n_sites, n_years];  // odds of detection
  
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all intervals
        //for(l in 1:n_visits){ // loop across all visits
        
          p[j,k] = inv_logit( // the inverse of the log odds of detection is equal to..
            p0 + # an intercept
            p_habitat_type * habitat_type[j] # a spatial detection effect
           ); // end p[j,k,l]
           
        //} // end loop across all visits
      } // end loop across all intervals
    } // end loop across all sites
  
  array[2, n_years] simplex[n_visits + 1] po; // Emission Probability
  
  
  for(j in 1:n_sites){
    for(k in 1 : n_years) {
      for(l in 1 : (n_visits + 1)) {
        po[1, k, l] = exp(binomial_lpmf(l - 1 | n_visits, p[j,k])); // something broke down here when trying to expand p
        po[2, k, l] = l == 1; 
      }
    }
  }

}

model {
  
  ////////////
  // PRIORS //
  ////////////
  
  // Flat priros Uniform(0, 1) are implicitly used on psi1, and ps.
  p0 ~ normal(0,2);
  p_habitat_type ~ normal(0,2);
  
  ////////////////
  // LIKELIHOOD //
  ////////////////
  
  // This implementation of the forward algorithm is derived from
  // Stan Modeling Language User's Guide and Reference Manual.
  
  for (j in 1 : n_sites) { // for each site
    array[2] real acc; // make an array of 2 numeric values 
    array[n_years] vector[2] gam; // and make an array of 2 numeric values for each year 
    
    // first year
    gam[1, 1] = psi1 * po[1, 1, sum_y[j, 1]]; // occupied * detections
    gam[1, 2] = (1 - psi1) * po[2, 1, sum_y[j, 1]]; // not occupied * detections
    
    // system dynamics
    for (k in 2 : n_years) {
      for (r in 1 : 2) {
        for (s in 1 : 2) {
          acc[s] = gam[k - 1, s] * ps[s, k - 1, r] * po[r, k, sum_y[j, k]];
        }
        gam[k, r] = sum(acc); // sum likelihood of state transitions for each year
      }
    }
    target += log(sum(gam[n_years])); // sum likelihood across all years
  }
  
} // end model

generated quantities {
  
  // Population occupancy, growth rate and turnover
  vector<lower=0, upper=1>[n_years] psi; // Occupancy probability
  vector<lower=0, upper=1>[ny_minus_1] phi; // Survival probability
  vector<lower=0, upper=1>[ny_minus_1] gamma; // Colonization probability
  array[n_sites, n_years] int<lower=0, upper=1> z; // Latent state of occurrence
  array[n_years] int<lower=0, upper=n_sites> n_occ; // Number of occupied sites
  vector[n_years - 1] growthr; // Population growth rate
  vector[n_years - 1] turnover; // Turnover rate
  
  // Latent state z[,] is estimated with a full simulation
  // unconditional on the observed V[,].
  for (k in 1 : ny_minus_1) {
    phi[k] = ps[1, k, 1];
    gamma[k] = ps[2, k, 1];
  }
  psi[1] = psi1; // psi in year one is quivalent to system initial state
  for (k in 2 : n_years) { // psi in following years is dynamically linked to initial state
    // psi[k-1] is a switch so if site was occupied in previous year, we estimate gamma
    // if empty in previous year we estimate colonization.
    psi[k] = psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1]; 
  }
  for (i in 1 : n_sites) {
    for (k in 1 : n_years) {
      z[i, k] = bernoulli_rng(psi[k]); // simulate z
    }
  }
  for (k in 1 : n_years) {
    n_occ[k] = sum(z[1 : n_sites, k]);
  }
  growthr[ : ny_minus_1] = psi[2 : ] ./ psi[ : ny_minus_1];
  turnover[ : ny_minus_1] = (1 - psi[ : ny_minus_1]) .* gamma[ : ny_minus_1]
                            ./ psi[2 : ];
}
