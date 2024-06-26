### Make figures from occupancy model outputs
# jcu; started oct 17, 2022

# library(rstan)
# library(bayesplot)
library(tidyverse)

stan_out <- readRDS("./analysis/stan_out.rds")
list_of_draws <- as.data.frame(stan_out)
print(stan_out, digits = 3)
colnames(list_of_draws)

n_species <- 32

fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,4] # 2.5
fit_summary$summary[1,8] # 97.5

#-------------------------------------------------
### Posterior model estimates for effect of habitat on detection and occupancy
X <- c(0.25, 0.75) # detection and abundance
Y <- c(mean(list_of_draws[,135]), mean(list_of_draws[,67])) # mean detection mean and mean abundance mean
# lower is mean - 1.96*sd
lower_95 <- c(fit_summary$summary[135,4], # detection
              fit_summary$summary[67,4])
upper_95 <- c(fit_summary$summary[135,8], # occupancy
              fit_summary$summary[67,8]) 

df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95))

species_estimates <- data.frame()

for(i in 1:n_species){
  
  species_estimates[1,i] <- mean(list_of_draws[,69+i]) # p
  species_estimates[2,i] <- mean(list_of_draws[,34+i]) # psi
  
}

df_estimates <- cbind(df_estimates, species_estimates)

(s <- ggplot(df_estimates) +
    geom_errorbar(aes(x=X, ymin=lower_95, ymax =upper_95),
                  color="black",width=0.1,size=1,alpha=0.5) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="Community Mean Effect of Meadow Habitat", breaks = c(0.25, 0.75),
                       labels=c("detection", "occurrence"),
                       limits = c(0,1)) +
    scale_y_continuous(name="Posterior Model Estimate",
                       limits = c(-3, 15)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

for(i in 1: n_species){
  
  test <- as.data.frame(cbind(X, df_estimates[1:2,4+i]))
  colnames(test) <- c("X", "Y")
  
  s <- s + geom_point(data = test, aes(x=X, y=Y), 
                      col = "skyblue", size = 12, shape = "-", alpha = 0.5)
  
}

(s <- s +
  geom_point(aes(x=X, y=Y),
             size = 6, alpha = 0.8) 
)

## --------------------------------------------------
### Effect of habitat on species occupancy
# plot showing habitat (x axis) and occupancy rate (y axis)
# with community mean trend line and uncertainty, 
# as well as species-specific trends

## --------------------------------------------------
## Prediction versus covariate
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

list_of_draws <- as.data.frame(stan_out)

# 95% conf int
plot(NA, xlim=c(-1,1), ylim=c(0,1),
     xlab = "Habitat enhancement",
     ylab = "Probability(Occurrence of the average species)")
curve(ilogit(fit_summary$summary[67,1]*x), add=TRUE, col = "blue", lwd = 3)
curve(ilogit(fit_summary$summary[67,4]*x), add=TRUE)
curve(ilogit(fit_summary$summary[67,8]*x), add=TRUE)

## --------------------------------------------------
## Natural habitat

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,67]
}

plot(NA, xlim=c(-1,1), ylim=c(0,1),
     xlab = "Habitat enhancement",
     ylab = "Pr(Occurrence)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[1,1] + # intercept
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[1,1] + # intercept
    fit_summary$summary[67,1]*x), 
  add=TRUE, col = "blue", lwd = 3)

low <- ilogit(
  fit_summary$summary[1,1] + # intercept
    (fit_summary$summary[67,1]*-2))

high <- ilogit(
  fit_summary$summary[1,1] + # intercept
    (fit_summary$summary[67,1]*2))

##---------------------------------------------------
### Prediction plot but with habitat enhancement as a categorical 
#-------------------------------------------------
### Posterior model estimates for effect of habitat on detection and occupancy
X <- c(0.25, 0.75) # occurrence control versus enhancement
Y <- c(
  # control
  ilogit(fit_summary$summary[1,1] + # intercept
    fit_summary$summary[67,1]*0),
  # enhanced
  ilogit(fit_summary$summary[1,1] + # intercept
           fit_summary$summary[67,1]*1)
  ) 

# lower is mean - 1.96*sd
lower_95 <- c(# control
  ilogit(fit_summary$summary[1,4] + # intercept
           fit_summary$summary[67,4]*0),
  # enhanced
  ilogit(fit_summary$summary[1,4] + # intercept
           fit_summary$summary[67,4]*1)
  )
upper_95 <- c(ilogit(fit_summary$summary[1,8] + # intercept
                       fit_summary$summary[67,8]*0),
              # enhanced
              ilogit(fit_summary$summary[1,8] + # intercept
                       fit_summary$summary[67,8]*1)
              ) 

df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95))

species_estimates <- data.frame()

for(i in 1:n_species){
  
  # control
  species_estimates[1,i] <- ilogit(mean(list_of_draws[,1+i]) + # species intercept
    mean(list_of_draws[,34+i])*0) # species slope (for ehancement = 0)
  # enhanced
  species_estimates[2,i] <- ilogit(mean(list_of_draws[,1+i]) + # species intercept
    mean(list_of_draws[,34+i])*1) # species slope (for ehancement = 1)
  
}

df_estimates <- cbind(df_estimates, species_estimates)

(p <- ggplot(df_estimates) +
    geom_errorbar(aes(x=X, ymin=lower_95, ymax =upper_95),
                  color="black",width=0.1,size=1,alpha=0.5) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="", breaks = c(0.25, 0.75),
                       labels=c("No enhancment", "Enhancement"),
                       limits = c(0,1)) +
    scale_y_continuous(name="Probability(Occurrence)",
                       limits = c(0, 1)) +
    guides(color = guide_legend(title = "")) +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

for(i in 1: n_species){
  
  test <- as.data.frame(cbind(X, df_estimates[1:2,4+i]))
  colnames(test) <- c("X", "Y")
  
  p <- p + geom_point(data = test, aes(x=X, y=Y), 
                      col = "skyblue", size = 24, shape = "-", alpha = 0.5)
  
}

(p <- p +
    geom_point(aes(x=X, y=Y),
               size = 6, alpha = 0.8) 
)



## --------------------------------------------------
### Effect of habitat on species detection
# plot showing habitat (x axis) and detection rate (y axis)
# with community mean trend line and uncertainty, 
# as well as species-specific trends

## --------------------------------------------------
### Species occupancy rate 
# geom_tile plot showing species specific occupancy intercepts
# ie how commonly does each species occur?