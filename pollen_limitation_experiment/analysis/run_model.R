## Analyze real data using a mixed-effects, logistic regression
#  started Dec. 14, 2022, J Ulrich

library(rstan)

# Format binary outcome data where.. 
# the outcome is: 
# 0 - a flower is not pollen limited, OR
# 1 - a flower is pollen limited
# A random intercept for planter pot, nested in site will account for heterogeneity
# in nested sample groups.

###-----------------------------------------------------------------------------
## Global options
max_PL_accepted = 2
PL_response_threshold = 0.5

source("./pollen_limitation_experiment/analysis/prep_data.R")
my_data <- prep_data(max_PL_accepted,
                     PL_response_threshold)


# data to feed to the model
N <- my_data$N # number of pairs
n_pots <- my_data$n_pots # number of pots
pot_ID <- as.numeric(as.factor(my_data$pots)) # pot ID names
n_sites <- my_data$n_sites # number of sites
site_chr <- my_data$sites
site_ID <- my_data$sites_numeric  # vector of sites
siteLookup <- as.numeric(as.factor(my_data$siteLookup)) 
herb_enhancement <- my_data$x # site type covariate
y <- my_data$y # outcome

herb_abundance_scaled <- my_data$herb_abundance_scaled

stan_data <- c("N", 
               # "n_pots", "pot_ID", "siteLookup",
               "n_sites", "site_ID",
               #"x", 
               "y", 
               #"herb_abundance_scaled"
               "herb_enhancement"
               )

# Parameters monitored
params <- c("beta0",
            "beta1",
            "beta0_site",
            "sigma_site",
            "sum_y_rep"
)

# MCMC settings
n_iterations <- 4000
n_thin <- 1
n_burnin <- 0.5*n_iterations
n_chains <- 4
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(beta0 = runif(1, -1, 1),
       sigma_site = runif(1, 0, 1),
       beta1 = runif(1, -1, 1)
       
  )
)

## --------------------------------------------------
### Run model

stan_model <- "./pollen_limitation_experiment/models/logistic_model_binary_covs.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores,
                     control = list(adapt_delta = 0.9))

print(stan_out, digits = 3)

saveRDS(stan_out, "./pollen_limitation_experiment/model_outputs/stan_out_logistic_regression.RDS")

stan_out <- readRDS("./pollen_limitation_experiment/model_outputs/stan_out_logistic_regression.RDS")


## --------------------------------------------------
### Simple diagnostic plots

print(stan_out, digits = 3)

# traceplot
traceplot(stan_out, pars = c(
  "beta0",
  "beta1",
  "sigma_site"
))

# pairs plot
pairs(stan_out, pars = c(
  "beta0",
  "beta1",
  "sigma_site"
))





## --------------------------------------------------
# posterior predictive check

fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


df_estimates <- data.frame(X = 1, 
                           Y = fit_summary$summary[15,6], 
                           lower_95 = fit_summary$summary[15,4],
                           upper_95 = fit_summary$summary[15,8],
                           lower_50 = fit_summary$summary[15,5],
                           upper_50 = fit_summary$summary[15,7]
) 




# posterior predictive check
# plot number of total pollen limited plants
df <- as.data.frame(cbind(herb_abundance_scaled, site_ID, y)) %>%
  mutate(pollen_limited_plants = sum(y)) %>%
  slice(1)

par(mar = c(4,5,5,2))
plot(1, type="n",
     xlim=c(0, 2), 
     xlab="",
     xaxt = "n",
     ylim=c(30,90), 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Real versus predicted number\n of pollen-limited plants")

# 95% quantiles
rect(xleft = (df_estimates$X-0.35), xright=(df_estimates$X+0.35), 
     ytop = df_estimates$lower_95, ybottom = df_estimates$upper_95,
     col = c_mid, border = NA
)

# 50% quantiles
rect(xleft =(df_estimates$X-0.35), xright=(df_estimates$X+0.35), 
     ytop = df_estimates$lower_50, ybottom = df_estimates$upper_50,
     col = c_mid_highlight, border = NA
)

# real observed value
points(x=df_estimates$X, y=df_estimates$Y, pch=19, cex=2)

## could also maybe do this by site if wanting more info
# group by site
# plot number of pollen limited plants
df <- as.data.frame(cbind(herb_abundance_scaled, site_ID, y)) %>%
  group_by(site_ID) %>%
  mutate(pollen_limited_plants = sum(y)) %>%
  ungroup() %>% 
  cbind(., site_chr) %>%
  group_by(site_ID) %>%
  slice(1) %>%
  ungroup()

# would need to write stan code to predict sums by site

# ...

## --------------------------------------------------
# prediction plot (continuous, see below for binary predictor)

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

tmp <- as.data.frame(stan_out) # take estimates from each HMC step as a df
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?
pred_length <- 100 # divide each covariate axis into some interval steps

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# get prediction sequence
original_herb_abundance <- my_data$original_herb_abundance
avg_log_flower_abundance <- my_data$avg_log_flower_abundance

original_herb <- seq(min(original_herb_abundance), 
                     max(original_herb_abundance), 
                     length.out = pred_length) # reasonable prediction range
herb_pred <- (original_herb - mean(original_herb_abundance)) / sd(original_herb_abundance) # unscale the dates


## --------------------------------------------------
# fill predictions for the expected outcome based on model param estimates

pred <- array(NA, dim=c(pred_length, n_samp)) # community means

for(i in 1:n_samp){
  # expected outcome for average site (no site random effect)
  pred[,i] <- ilogit( 
    # herbaceous flower trend
    # beta0 +
    tmp[i,1] + 
    # beta1*herb_pred +
    tmp[i,2]*herb_pred
  )
}

# posterior means by community average 
cri <- apply(pred, 
             c(1), # apply across first dimension (for each herb value, find the mean and CI)
             function(x) quantile(x, 
             prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

herb_df <- as.data.frame(cbind(original_herb, cri[3,], 
                               cri[1,], cri[5,],
                               cri[2,], cri[4,])) %>%
  rename("mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6")

# add observed values
tmp2 <- as.data.frame(cbind(avg_log_flower_abundance, y))

p <- ggplot(data = herb_df, aes(original_herb, mean)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=lower95, 
    ymax=upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  geom_jitter(data = tmp2, aes(x=avg_log_flower_abundance, y=y),
              width = 0, height = 0.05, alpha=0.5, size=2) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(-0.05, 1.05)) +
  theme_bw() +
  ylab("Pr(pollen limitation)") +
  xlab("log(herbaceous flower abundance / 20m^2)") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

## --------------------------------------------------
# prediction plot (binary)

fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

tmp <- as.data.frame(stan_out) # take estimates from each HMC step as a df
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?
pred_length <- 2 # divide each covariate axis into some interval steps

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# get prediction sequence
original_herb <- c(0,1)
herb_pred <- c(0,1)

## --------------------------------------------------
# fill predictions for the expected outcome based on model param estimates

pred <- array(NA, dim=c(pred_length, n_samp)) # community means

for(i in 1:n_samp){
  # expected outcome for average site (no site random effect)
  pred[,i] <- ilogit( 
    # herbaceous flower trend
    # beta0 +
    tmp[i,1] + 
      # beta1*herb_pred +
      tmp[i,2]*herb_pred
  )
}

# posterior means by community average 
cri <- apply(pred, 
             c(1), # apply across first dimension (for each herb value, find the mean and CI)
             function(x) quantile(x, 
                                  prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

herb_df <- as.data.frame(cbind(original_herb, cri[3,], 
                               cri[1,], cri[5,],
                               cri[2,], cri[4,])) %>%
  rename("mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6")

# add observed values
tmp2 <- as.data.frame(cbind(herb_enhancement, y))

p <- ggplot(data = herb_df, aes(original_herb, mean)) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=lower95, 
    ymax=upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  geom_jitter(data = tmp2, aes(x=herb_enhancement, y=y),
              width = 0.1, height = 0.1, alpha=0.5, size=2.5) +
  xlim(c(min(original_herb), max(original_herb))) +
  ylim(c(-0.05, 1.05)) +
  theme_bw() +
  ylab("Pr(pollen limitation)") +
  xlab("") +
  scale_x_continuous(limits = c(-0.15, 1.15),
                     breaks = c(0, 1),
                     labels = c("control","herb. enhancement")) +
  scale_y_continuous(limits = c(-0.1, 1.1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p


## --------------------------------------------------
# could also do this with canned rstanarm package
fit1 <- rstanarm::stan_glmer(y ~ herb_abundance_scaled + (1|site_ID), 
                     family = binomial(link = "logit"))

rstanarm::pp_check(fit1)


hist(y)
