// dynamic multi-species occupancy model for pollinators in urban parks
// jcu, started nov 23, 2023.
// see work done here for some ideas on translating JAGS models into STAN:
// https://github.com/stan-dev/example-models/blob/master/BPA/Ch.13/Dynocc.stan#L74
// there are some tricky things M. Joseph has done to build the model template
// will need to walk through to understand fully

data {
  
  int<lower=1> n_sites;  // sites within region
  //int<lower=1> sites[n_sites];  // vector of sites
  
  int<lower=1> n_years; // number of years surveyed
  int<lower=1> n_visits; // number of MAX surveys per year
  
  int<lower=0> V[n_sites, n_years, n_visits]; // visits l when species i was detected at site j on interval k
  
} // end data

transformed data {
  array[n_sites, n_years] int<lower=1, upper=n_visits + 1> sum_y; // sum + 1
  // calculate total years with transitions
  int ny_minus_1 = n_years - 1;
  
  for (j in 1 : n_sites) {
    for (k in 1 : n_years) {
      sum_y[j, k] = sum(V[j, k, 1 : n_visits]) + 1;
    }
  }
} // end  transformed data

parameters {
  real<lower=0, upper=1> psi1; // Occupancy probability at t=1
  vector<lower=0, upper=1>[n_years] p; // Detection probability
  array[2, n_years - 1] simplex[2] ps; // Transition probability
  // This is equivalent to the following.
  //  ps[1, t, 1] = phi[t];
  //  ps[1, t, 2] = 1.0 - phi[t];
  //  ps[2, t, 1] = gamma[t];
  //  ps[2, t, 2] = 1.0 - gamma[t];
}

transformed parameters {
  array[2, n_years] simplex[n_visits + 1] po; // Emission Probability
  
  for (t in 1 : n_years) {
    for (r in 1 : (n_visits + 1)) {
      po[1, t, r] = exp(binomial_lpmf(r - 1 | n_visits, p[t])); // occupied
      po[2, t, r] = r == 1; // not occupied
    }
  }
}

model {
  
  ////////////
  // PRIORS //
  ////////////
  
  // Flat priros Uniform(0, 1) are implicitly used on psi1, p and ps.
  
  ////////////////
  // LIKELIHOOD //
  ////////////////
  
  // This implementation of the forward algorithm is derived from
  // Stan Modeling Language User's Guide and Reference Manual.
  for (i in 1 : n_sites) { // for each site
    array[2] real acc; // make an array of...
    array[n_years] vector[2] gam;
    
    // First year
    gam[1, 1] = psi1 * po[1, 1, sum_y[i, 1]];
    gam[1, 2] = (1 - psi1) * po[2, 1, sum_y[i, 1]];
    
    for (t in 2 : n_years) {
      for (k in 1 : 2) {
        for (j in 1 : 2) {
          acc[j] = gam[t - 1, j] * ps[j, t - 1, k] * po[k, t, sum_y[i, t]];
        }
        gam[t, k] = sum(acc);
      }
    }
    target += log(sum(gam[n_years]));
  }
  
} // end model

generated quantities {
  // Population occupancy, growth rate and turnover
  vector<lower=0, upper=1>[n_years] psi; // Occupancy probability
  vector<lower=0, upper=1>[ny_minus_1] phi; // Survival probability
  vector<lower=0, upper=1>[ny_minus_1] gamma; // Colonization probability
  array[n_sites, n_years] int<lower=0, upper=1> z; // Latent state of occurrence
  array[n_years] int<lower=0, upper=n_sites> n_occ; // Number of occupancy
  vector[n_years - 1] growthr; // Population growth rate
  vector[n_years - 1] turnover; // Turnover rate
  
  // Latent state z[,] is estimated with a full simulation
  // unconditional on the observed V[,].
  for (k in 1 : ny_minus_1) {
    phi[k] = ps[1, k, 1];
    gamma[k] = ps[2, k, 1];
  }
  psi[1] = psi1;
  for (k in 2 : n_years) {
    psi[k] = psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1];
  }
  for (i in 1 : n_sites) {
    for (k in 1 : n_years) {
      z[i, k] = bernoulli_rng(psi[k]);
    }
  }
  for (t in 1 : n_years) {
    n_occ[t] = sum(z[1 : n_sites, t]);
  }
  growthr[ : ny_minus_1] = psi[2 : ] ./ psi[ : ny_minus_1];
  turnover[ : ny_minus_1] = (1 - psi[ : ny_minus_1]) .* gamma[ : ny_minus_1]
                            ./ psi[2 : ];
}
