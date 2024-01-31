// dynamic multi-species occupancy model for pollinators in urban parks
// jcu, started nov 23, 2023.
// see work done here for some ideas on translating JAGS models into STAN:
// https://www.r-bloggers.com/2014/11/dynamic-occupancy-models-in-stan/

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

// multinormal model with correlated species effects for all 4 levels, with 6 rho's

functions {
  
  // covariance matrix for detection rates from different data sources
  matrix custom_cov_matrix(vector sigma, real rho1, real rho2, real rho3,
                            real rho4, real rho5, real rho6) {
                              
    matrix[4,4] Sigma;
    
    Sigma[1,1] = square(sigma[1]); // species variation in occurrence1 rates
    Sigma[2,2] = square(sigma[2]); // species variation in gamma rates
    Sigma[3,3] = square(sigma[2]); // species variation in phi rates
    Sigma[4,4] = square(sigma[2]); // species variation in detection rates
    
    Sigma[1,2] = sigma[1] * sigma[2] * rho1; // correlation between species-specific occurrence1 and gamma rates
    Sigma[2,1] = Sigma[1,2];    
    
    Sigma[1,3] = sigma[1] * sigma[3] * rho2; // correlation between species-specific occurrence1 and phi rates
    Sigma[3,1] = Sigma[1,3];    
    
    Sigma[1,4] = sigma[1] * sigma[4] * rho3; // correlation between species-specific occurrence1 and detection rates
    Sigma[4,1] = Sigma[1,4]; 
    
    Sigma[2,3] = sigma[2] * sigma[3] * rho4; // correlation between species-specific gamma and phi rates
    Sigma[3,2] = Sigma[2,3]; 
    
    Sigma[2,4] = sigma[2] * sigma[4] * rho5; // correlation between species-specific gamma and detection rates
    Sigma[4,2] = Sigma[2,4]; 
    
    Sigma[3,4] = sigma[3] * sigma[4] * rho6; // correlation between species-specific phi and detection rates
    Sigma[4,3] = Sigma[3,4]; 
    
    return Sigma;
  }
  
  // mean values for multivariate normal distribution (center species on global intercepts)
  vector mu(real psi1_0, real gamma0, real phi0, real p0){
    vector[4] global_intercepts;
    global_intercepts[1] = psi1_0;
    global_intercepts[2] = gamma0;
    global_intercepts[3] = phi0;
    global_intercepts[4] = p0;
    return global_intercepts;
  }
  
}

data {
  
  int<lower=0> n_species;
  int<lower=1> species[n_species]; // vector of species identities
  int<lower=0> n_sites;
  int<lower=1> sites[n_sites]; // vector of site identities
  int<lower=0> n_years;
  int<lower=0> n_visits;
  int<lower=0,upper=1> V[n_species, n_sites, n_years, n_visits];
   
  // covariate data
  int<lower=0, upper=1> habitat_type[n_sites]; // categorical habitat type (0 or 1)
  real date_scaled[n_sites, n_years, n_visits]; // scaled day of year on which each visit was conducted

} // end data

parameters {
  
  // Covararying Parameters
  real<lower=-1,upper=1> rho1;  // correlation of psi1 gamma 
  real<lower=-1,upper=1> rho2;  // correlation of... psi1 phi
  real<lower=-1,upper=1> rho3;  // correlation of... psi1 detection
  real<lower=-1,upper=1> rho4;  // correlation of... gamma phi
  real<lower=-1,upper=1> rho5;  // correlation of... gamma detecion
  real<lower=-1,upper=1> rho6;  // correlation of... phi detection
  vector<lower=0>[4] sigma_species; // variance in species-specific detection rates 
  vector[4] species_intercepts[n_species];// species-level detection intercepts
  
  real psi1_0;
  //vector[n_species] psi1_species;
  //real sigma_psi1_species;
  vector[n_species] psi1_habitat;
  real mu_psi1_habitat;
  real sigma_psi1_habitat;
  real gamma0;
  //vector[n_species] gamma_species;
  //real sigma_gamma_species;
  vector[n_species] gamma_habitat;
  real mu_gamma_habitat;
  real sigma_gamma_habitat;
  real phi0;
  //vector[n_species] phi_species;
  //real sigma_phi_species;
  vector[n_species] phi_habitat;
  real mu_phi_habitat;
  real sigma_phi_habitat;
  
  real p0;
  //vector[n_species] p_species;
  //real sigma_p_species;
  //vector[n_sites] p_site;
  //real sigma_p_site;
  real p_habitat;
  vector[n_species] p_date;
  real mu_p_species_date;
  real sigma_p_species_date;
  vector[n_species] p_date_sq;
  real sigma_p_species_date_sq;
  real mu_p_species_date_sq;

} // end parameters

transformed parameters {
   
  // logit scale psi1, gamma, phi
  real psi1[n_species, n_sites]; // odds of occurrence year 1
  real gamma[n_species, n_sites, n_years]; // odds of colonization
  real phi[n_species, n_sites, n_years]; // odds of persistence
   
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          //psi1_species[species[i]] + // a species specific intercept
          species_intercepts[species[i],1] + // species random effect
          psi1_habitat[species[i]] * habitat_type[j] // a spatial effect
          ); // end phi[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          //gamma_species[species[i]] + // a species specific intercept
          species_intercepts[species[i],2] + // species random effect
          gamma_habitat[species[i]] * habitat_type[j] // a spatial effect
          ); // end phi[j,k]
        
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          //phi_species[species[i]] + // a species specific intercept
          species_intercepts[species[i],3] + // species random effect
          phi_habitat[species[i]] * habitat_type[j] // a spatial effect
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
            //p_species[species[i]] + // a species-specific intercept
            species_intercepts[species[i],4] + // species random effect
            //p_site[sites[j]] + // a site-specific intercept
            p_habitat * habitat_type[j] + // a spatial detection effect
            p_date[species[i]] * date_scaled[j,k,l] + // a species-specific phenological detection effect (peak)
            p_date_sq[species[i]] * (date_scaled[j,k,l])^2 // a species-specific phenological detection effect (decay)
            ); // end p[j,k,l]
             
        } // end loop across all visits
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
  
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // correlated params:
  // correlated species effects
  sigma_species[1] ~ normal(0, 2); // initial occurrence
  sigma_species[2] ~ normal(0, 1); // colonization
  sigma_species[3] ~ normal(0, 1); // persistence
  sigma_species[4] ~ normal(0, 2); // detection
  (rho1 + 1) / 2 ~ beta(2, 2);
  (rho2 + 1) / 2 ~ beta(2, 2);
  (rho3 + 1) / 2 ~ beta(2, 2);
  (rho4 + 1) / 2 ~ beta(2, 2);
  (rho5 + 1) / 2 ~ beta(2, 2);
  (rho6 + 1) / 2 ~ beta(2, 2);
  
  // correlated species-specific detection rates
  // will send the mean (mu), variance and correlation to the covariance matrix
  species_intercepts ~ multi_normal(mu(psi1_0, gamma0, phi0, p0), 
    custom_cov_matrix(sigma_species, rho1, rho2, rho3, rho4, rho5, rho6));
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0,1); // initial occupancy rate
  //psi1_species ~ normal(psi1_0, sigma_psi1_species); // species-specific intercepts (centered on global)
  //sigma_psi1_species ~ normal(0, 1); // variation in species-specific intercepts
  psi1_habitat ~ normal(mu_psi1_habitat,sigma_psi1_habitat); // effect of habitat on occurrence
  mu_psi1_habitat ~ normal(0,2); // community mean
  sigma_psi1_habitat ~ normal(0,1); // species variation
  // colonization
  gamma0 ~ normal(0,1); // colonization rate
  //gamma_species ~ normal(gamma0, sigma_gamma_species); // species-specific intercepts (centered on global)
  //sigma_gamma_species ~ normal(0, 1); // variation in species-specific intercepts
  gamma_habitat ~ normal(mu_gamma_habitat,sigma_gamma_habitat); // effect of habitat on colonization
  mu_gamma_habitat ~ normal(0,2); // community mean
  sigma_gamma_habitat ~ normal(0,1); // species variation
  // persistence
  phi0 ~ normal(0,1); // global persistence intercept
  //phi_species ~ normal(phi0, sigma_phi_species); // species-specific intercepts (centered on global)
  //sigma_phi_species ~ normal(0, 1); // variation in species-specific intercepts
  phi_habitat ~ normal(mu_phi_habitat, sigma_phi_habitat); // effect of habitat on persistence
  mu_phi_habitat ~ normal(0,2); // community mean
  sigma_phi_habitat ~ normal(0,1); // species variation
  
  // detection
  p0 ~ normal(0,1); // global intercept
  //p_species ~ normal(p0, sigma_p_species); // species-specific intercepts (centered on global)
  //sigma_p_species ~ normal(0, 1); // variation in species-specific intercepts
  //p_site ~ normal(0, sigma_p_site); // site-specific intercepts
  //sigma_p_site ~ normal(0,1); // variation in site-specific intercepts
  p_habitat ~ normal(0,2); // effect of habitat on detection
  // phenology X detection
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0,2); // mean
  sigma_p_species_date ~ normal(0, 1); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0,2); // mean
  sigma_p_species_date_sq ~ normal(0, 1); // variation
  
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
  
  // Diversity estimation
  // number of species at each site in each year
  int z_simmed[n_species, n_sites, n_years]; // simulate occurrence
  int species_richness[n_sites, n_years]; // site/year species richness
  real avg_species_richness_enhanced[n_years]; // average across sites
  real avg_species_richness_control[n_years]; // average across sites
  
  // equilibrium occupancy rate (for average species)  
  real psi_eq_habitat0; // control sites
  real psi_eq_habitat1; // enhanced sites
  
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
  
  psi_eq_habitat0 = inv_logit(gamma0) / (inv_logit(gamma0)+(1-inv_logit(phi0))); // equilibrium occupancy rate 
  psi_eq_habitat1 = inv_logit(gamma0 + mu_gamma_habitat) / (inv_logit(gamma0 + mu_gamma_habitat)+(1-inv_logit(phi0 + mu_phi_habitat))); // equilibrium occupancy rate
  
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
            
            Z[i,j,k] = bernoulli_logit_rng(psi[i,j,k]);
            
            // occupancy but never observed by either dataset
            //real ulo = inv_logit(psi[i,j,k]) * 
              // ((1 - inv_logit(p[j,k]))^n_visits) // for no visit heterogeneity in p
              //(1 - inv_logit(p[i,j,k,1])) * (1 - inv_logit(p[i,j,k,2])) *
              //(1 - inv_logit(p[i,j,k,3])) * (1 - inv_logit(p[i,j,k,4])) *
              //(1 - inv_logit(p[i,j,k,5])) * (1 - inv_logit(p[i,j,k,6]));
            // non-occupancy
           //real uln = (1 - inv_logit(psi[i,j,k]));
            
            // outcome of occupancy given the likelihood associated with both possibilities
            //Z[i,j,k] = bernoulli_rng(ulo / (ulo + uln));
            
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
