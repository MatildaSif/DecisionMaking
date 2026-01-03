pacman::p_load(extraDistr, R2jags, parallel, ggpubr, truncnorm)

setwd("/work/JoMat/DecisionMaking/src/VSE")

print("About to source heir_VSE_sim.R...")
flush.console()

source('heir_VSE_sim.R')

print("Successfully loaded heir_VSE_sim.R!")
flush.console()

# Try running it once with dummy data
print("Testing heir_VSE_sim function...")
flush.console()

R_test <- matrix(runif(48*100, 0, 100), nrow=48, ncol=100)
L_test <- matrix(runif(48*100, 0, 100), nrow=48, ncol=100)

result <- heir_VSE_sim(R_test, L_test, nsubs=2, ntrials_all=c(100, 100),
                       mu_theta=0.5, mu_delta=0.7, mu_alpha=0.3, 
                       mu_phi=0.4, mu_c=0.8,
                       sigma_theta=0.1, sigma_delta=0.1, 
                       sigma_alpha=0.1, sigma_phi=0.1, sigma_c=0.1)

print("Successfully ran heir_VSE_sim!")
flush.console()