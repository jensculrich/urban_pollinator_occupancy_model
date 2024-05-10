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
  
  vector[N] x; // covariate of site type for each flower pair
  
  int<lower=0,upper=1> y[N]; // Pollen limitation (PL) outcome for each flower pair
  
  int<lower=1> n_pots;  // number of pots // number of level-2 clusters
  int<lower=1, upper=n_pots> pot_ID[N];  // vector of pot names // level-2 clusters
  
  int<lower=1> n_sites;  // number of sites // number of level-3 clusters
  int<lower=1, upper=n_sites> site_ID[N];  // vector of site names // level-3 clusters
  
  int<lower=1> siteLookup[n_pots]; // level-3 cluster look up vector for level-2 cluster
}

parameters {
  
  real alpha0; // global intercept 
  
  // Level-2 random effect
  // pot specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all pots.
  real u_pot[n_pots]; // pot specific intercept for PL outcome
  real<lower=0> sigma_alpha_pot; // variance in pot intercepts
  
  // Level-3 random effect
  // site specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all sites.
  real u_site[n_sites]; // site specific intercept for PL outcome
  real<lower=0> sigma_alpha_site; // variance in site intercepts
  
  real beta; // effect of site type on the PL outcome

}

transformed parameters{
  
  // varying intercepts
  real alpha0_pot[n_pots];
  real alpha0_site[n_sites];

  // the linear predictor for the individual observations
  real p[N];

  // compute the varying intercept at the site level
  // Level-3 (n_sites level-3 random intercepts)
  for(i in 1:n_sites){
    alpha0_site[i] = u_site[i];
  }

  // compute varying intercept at the pot within site level
  // Level-2 (n_pots level-2 random intercepts)
  for(i in 1:n_pots){
     alpha0_pot[i] = alpha0_site[siteLookup[i]] + u_pot[i];
  }

  // Individual flower mean
  for(i in 1:N){
    
      p[i] = alpha0 + // a global intercept
             alpha0_pot[pot_ID[i]] + // a pot|site specific intercept
             beta * x[i] // an effect of site type
            ; // end p
              
  }
  
}

model {
  
  // PRIORS
  
  alpha0 ~ cauchy(0, 2.5); // weakly informative prior for global intercept
  
  
  // level-3 grouping
  u_site ~ normal(0, sigma_alpha_site); 
  // prob of success intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_alpha_site ~ cauchy(0, 1); // weakly informative prior
  
  // level-2 grouping
  u_pot ~ normal(0, sigma_alpha_pot); 
  // prob of success intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_alpha_pot ~ cauchy(0, 1); // weakly informative prior
  
  beta ~ cauchy(0, 2.5); // weakly informative prior for effect of site type on outcome
  
  // LIKELIHOOD
  
  y ~ bernoulli_logit(p);
  
}
