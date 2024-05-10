library(rstanarm)

stan_out <- readRDS("./model_outputs/stan_out_logistic_regression.RDS")

# To do the ppc with this model output, should write a data prediction into 
# the generated quantities block. Could do it manually using the fit matrix but 
# it could be a lot more work and easier to have it output with the model.
fit_matrix <- as.matrix(stan_out)

###------------------------------------------------------------------------------
## Posterior Predictive Check

# can the model can reproduce the mean observed in the real data?
# (where the mean is the proportion of plants that are 
# pollen limited accourding to our pollen limitation threshold)
(test <- pp_check(stan_out, plotfun = "stat", binwidth = 0.01))

# The histogram shows the distributions of the fraction of simulated data 
# that resulted in pollen limitation in each of the posterior draws.

# can the model can reproduce the patterns observed in the real data
# where we have most plants being pollinated relatively sufficiently
# and about 1/3 being pollen limited?
yrep <- posterior_predict(stan_out, draws = 500)
pp_check(stan_out, nreps=100)