	
data {
   int<lower=0> n_sites;
   int<lower=0> n_years;
   int<lower=0> n_visits;
   int<lower=0,upper=1> V[n_sites, n_years, n_visits];
}
parameters {
   real<lower=0,upper=1> p;
   real<lower=0,upper=1> gamma;
   real<lower=0,upper=1> phi;
   real<lower=0, upper=1> psi1;
}
transformed parameters {
   matrix[n_sites, n_years] psi;
   for (j in 1:n_sites){
     for (k in 1:n_years){
       if (k < 2){ // initial state
          psi[j, k] = psi1; 
       } else { // system dynamics
          psi[j, k] = psi[j, k-1] * phi + (1 - psi[j, k-1]) * gamma; 
       }
     }
   }
}
model {
   // priors
  psi1 ~ uniform(0,1);
  gamma ~ uniform(0,1);
  phi ~ uniform(0,1);
  p ~ uniform(0,1);
  
   // likelihood
  for (j in 1:n_sites){
    for (k in 1:n_years){
      if (sum(V[j, k]) > 0){
        target += (log(psi[j, k]) + bernoulli_lpmf(V[j, k]|p));
      } else {
        target += (log_sum_exp(log(psi[j, k]) + bernoulli_lpmf(V[j, k]|p),
                                      log(1-psi[j, k])));
      }
    }
  }
}
