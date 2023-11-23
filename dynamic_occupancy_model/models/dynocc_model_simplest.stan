	
data {
   int<lower=0> nsite;
   int<lower=0> nyear;
   int<lower=0> nrep;
   int<lower=0,upper=1> Y[nsite,nyear,nrep];
}
parameters {
   real<lower=0,upper=1> p;
   real<lower=0,upper=1> gamma;
   real<lower=0,upper=1> phi;
   real<lower=0, upper=1> psi1;
}
transformed parameters {
   matrix[nsite, nyear] psi;
   for (r in 1:nsite){
     for (t in 1:nyear){
       if (t < 2){
          psi[r, t] <- psi1;
       } else {
          psi[r, t] <- psi[r, t-1] * phi + (1 - psi[r, t-1]) * gamma;
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
  for (r in 1:nsite){
    for (t in 1:nyear){
      if (sum(Y[r, t]) > 0){
        increment_log_prob(log(psi[r, t]) + bernoulli_log(Y[r, t], p));
      } else {
        increment_log_prob(log_sum_exp(log(psi[r, t]) + bernoulli_log(Y[r, t],p),
                                      log(1-psi[r, t])));
      }
    }
  }
}
