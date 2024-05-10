// STAN MODEL FOR LOGISTIC MIXED-EFFECTS REGRESSION
// started Dec. 14, 2022, J Ulrich

// Model is intended to estimate the effect of categorical site type
// on the outcome of: 
// 0 - a flower is not pollen limited, OR
// 1 - a flower is pollen limited
// A random intercept for plants grouped in planter pots, 
// nested in sites, will account for heterogeneity among sample groups.

data {
  
  int<lower=0> N; // number of paired flowers
  
  int<lower=0,upper=1> y[N]; // Pollen limitation (PL) outcome for each flower pair
  
  int<lower=1> n_sites;  // number of sites // number of level-3 clusters
  int<lower=1, upper=n_sites> site_ID[N];  // vector of site names // level-3 clusters
  
  vector[N] herb_abundance_scaled; // continuous herb flower abundance predictor
  
}

parameters {
  
  real beta0; // global intercept 
  
  // Level-3 random effect
  // site specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all sites.
  real beta0_site[n_sites]; // site specific intercept for PL outcome
  real<lower=0> sigma_site; // variance in site intercepts
  
  real beta1; // effect of site type on the PL outcome

}

transformed parameters{

  // the linear predictor for the individual observations
  // prob. of being pollen limited
  real p[N];

  // Individual flower mean
  for(i in 1:N){
    
      p[i] = beta0 + // a global intercept
             beta0_site[site_ID[i]] + // a pot|site specific intercept
             beta1 * herb_abundance_scaled[i] // an effect of site type
            ; // end p
              
  }
  
}

model {
  
  // PRIORS
  
  beta0 ~ normal(0, 2); // weakly informative prior for global intercept
  
  
  // level-3 grouping
  beta0_site ~ normal(0, sigma_site); 
  // prob of success intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_site ~ cauchy(0, 1); // weakly informative prior

  beta1 ~ normal(0, 2); // weakly informative prior for effect of site type on outcome
  
  // LIKELIHOOD
  
  y ~ bernoulli_logit(p);
  
}

generated quantities {
  
  // posterior predictive check
  // group by site and predict the number of pollen limited plants
  // compare the number of pollen limited plants at each site to the
  // number of pollen limited plants in the real data
  
  int<lower=0> y_rep[N]; // simulated PL outcomes
  int<lower=0> sum_y_rep; // total pollen limited plants
  
  // generating posterior predictive distribution
  // Predict pollen limitation at each plant,
  for(i in 1:N) { // loop across all plants
    y_rep[i] = bernoulli_logit_rng(p[i]);
  } 
  
  sum_y_rep = sum(y_rep);
  
}
