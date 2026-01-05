install.packages("pacman")
#pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)
install.packages("ggpubr")
install.packages("extraDistr")
install.packages("truncnorm")
install.packages("R2jags")
library(ggpubr)
library(R2jags)
library(parallel)
library(extraDistr)
library(truncnorm)

set.seed(1983)

setwd('/work/JoMat/DecisionMaking/src/PVL/')

# Create output folder structure
output_dir <- "../../PVL/outputs/parameter_recovery"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 10 # size of our subdivisions for pseudorandomization
freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
bad_r <- 100 # "bad" winnings
bad_freq_l <- -250 # "bad" frequent loss
bad_infreq_l <- -1250 # "bad" infrequent loss
good_r <- 50 # "good" winnings
good_freq_l <- -50 # "good" frequent loss
good_infreq_l <- -250 # "good" infrequent loss

# Bad frequent
A_R <- rep(bad_r, nstruct) # we win on every trials
A_L <- c(rep(bad_freq_l, nstruct*freq),rep(0,nstruct*(1-freq))) # we have losses half of the time

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq))) # we have losses 1/10th of the time

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq),rep(0,nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq)))

# create the pseudorandomized full payoff structure
A <- array(NA,ntrials) # setting up and empty array to be filled
B <- array(NA,ntrials)
C <- array(NA,ntrials)
D <- array(NA,ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}

# "less generic" code
# A <- c()
# B <- c()
# C <- c()
# D <- c()
# for (i in 1:10) {
#   A <- (append(A,A_R + sample(A_L)))
#   B <- (append(B,B_R + sample(B_L)))
#   C <- (append(C,C_R + sample(C_L)))
#   D <- (append(D,D_R + sample(D_L)))
# }

payoff <- cbind(A,B,C,D)/100 # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)

###--------------Run full parameter recovery -------------
niterations <- 10 # was 100 # fewer because it takes too long
nsubs <- 48 # mimicking the data structure from Ahn et al
ntrials_all <- rep(100, 48) # all 48 subs have 100 trials each

# mu
true_mu_w <- array(NA,c(niterations))
true_mu_A <- array(NA,c(niterations))
true_mu_theta <- array(NA,c(niterations))
true_mu_a <- array(NA,c(niterations))

infer_mu_w <- array(NA,c(niterations))
infer_mu_A <- array(NA,c(niterations))
infer_mu_theta <- array(NA,c(niterations))
infer_mu_a <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_sigma_w <- array(NA,c(niterations))
true_sigma_A <- array(NA,c(niterations))
true_sigma_theta <- array(NA,c(niterations))
true_sigma_a <- array(NA,c(niterations))

infer_lambda_w <- array(NA,c(niterations))
infer_lambda_A <- array(NA,c(niterations))
infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_a <- array(NA,c(niterations))

start_time = Sys.time()
for (i in 1:niterations) {
  ntrials <- ntrials_all
  
  mu_w <- runif(1,.5,2.5)
  mu_A <- runif(1,0,1)
  mu_theta <- runif(1, 0.5, 1.5)
  mu_a <- runif(1,0,1)
  
  sigma_w <- runif(1,0.01,0.2) # set all minimums to 0.01 to avoid dividing by 0 --> get infinities
  sigma_A <- runif(1,0.01,0.1)
  sigma_theta <- runif(1, 0.05, 0.15)
  sigma_a <- runif(1,0.01,0.1)
  
  source('hier_PVL_sim.R')
  PVL_sims <- hier_PVL_sim(payoff,nsubs,ntrials,mu_w,mu_A,mu_a,mu_theta,
                  sigma_w, sigma_A, sigma_a, sigma_theta)
  
  x <- PVL_sims$x
  X <- PVL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","ntrials","nsubs") 
  params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="hier_PVL.txt", n.chains=3, 
                           n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # mu
  true_mu_w[i] <- mu_w
  true_mu_A[i] <- mu_A
  true_mu_theta[i] <- mu_theta
  true_mu_a[i] <- mu_a
  
  # sigma
  true_sigma_w[i] <- sigma_w
  true_sigma_A[i] <- sigma_A
  true_sigma_theta[i] <- sigma_theta
  true_sigma_a[i] <- sigma_a
  

  # posterior mean (simplest and robust)
  Y <- samples$BUGSoutput$sims.list
  infer_mu_w[i] <- mean(Y$mu_w)
  infer_mu_A[i] <- mean(Y$mu_A)
  infer_mu_theta[i] <- mean(Y$mu_theta)
  infer_mu_a[i] <- mean(Y$mu_a)
  
  # similarly for lambda / sigma
  infer_lambda_w[i] <- mean(Y$lambda_w)
  infer_lambda_A[i] <- mean(Y$lambda_A)
  infer_lambda_theta[i] <- mean(Y$lambda_theta)
  infer_lambda_a[i] <- mean(Y$lambda_a)
  
  print(i)
  
}

end_time = Sys.time()
end_time - start_time

# scatter plots
source('recov_plot.R')
pl1 <- recov_plot(true_mu_w, infer_mu_w, c("true mu_w", "infer mu_w"), 'smoothed linear fit')
pl2 <- recov_plot(true_mu_A, infer_mu_A, c("true mu_A", "infer mu_A"), 'smoothed linear fit')
pl3 <- recov_plot(true_mu_theta, infer_mu_theta, c("true mu_theta", "infer mu_theta"), 'smoothed linear fit')
pl4 <- recov_plot(true_mu_a, infer_mu_a, c("true mu_a", "infer mu_a"), 'smoothed linear fit')

# Arrange and display
mu_plots <- ggarrange(pl1, pl2, pl3, pl4, ncol = 2, nrow = 2)
print(mu_plots)

# Save mu parameters plot
ggsave(file.path(output_dir, "PVL_recovery_mu_parameters.png"), mu_plots, 
       width = 12, height = 10, dpi = 300)
cat("\nSaved: PVL_recovery_mu_parameters.png\n")

# convert sigma to precision
true_lambda_w <- 1 / (true_sigma_w^2)
true_lambda_A <- 1 / (true_sigma_A^2)
true_lambda_theta <- 1 / (true_sigma_theta^2)
true_lambda_a <- 1 / (true_sigma_a^2)

# sigma (aka. true_lambda) re-coded as precision
pl5 <- recov_plot(infer_lambda_w, 1/(true_lambda_w^2), c("infer lambda_w","true lambda_w"), 'smoothed linear fit')
pl6 <- recov_plot(infer_lambda_A, 1/(true_lambda_A^2), c("infer lambda_A","true lambda_A"), 'smoothed linear fit')
pl7 <- recov_plot(infer_lambda_theta, 1/(true_lambda_theta^2), c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl8 <- recov_plot(infer_lambda_a, 1/(true_lambda_a^2), c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')

# Arrange and display
lambda_precision_plots <- ggarrange(pl5, pl6, pl7, pl8, ncol = 2, nrow = 2)
print(lambda_precision_plots)

# Save lambda precision plot
ggsave(file.path(output_dir, "PVL_recovery_lambda_precision.png"), lambda_precision_plots, 
       width = 12, height = 10, dpi = 300)
cat("Saved: PVL_recovery_lambda_precision.png\n")


# lambda (aka. infer_lambda) re-coded as SD
pl9 <- recov_plot(1/sqrt(infer_lambda_w), true_lambda_w, c("infer lambda_w","true lambda_w"), 'smoothed linear fit')
pl10 <- recov_plot(1/sqrt(infer_lambda_A), true_lambda_A, c("infer lambda_A","true lambda_A"), 'smoothed linear fit')
pl11 <- recov_plot(1/sqrt(infer_lambda_theta), true_lambda_theta, c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl12 <- recov_plot(1/sqrt(infer_lambda_a), true_lambda_a, c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')

# Arrange and display
lambda_sd_plots <- ggarrange(pl9, pl10, pl11, pl12, ncol = 2, nrow = 2)
print(lambda_sd_plots)

# Save lambda SD plot
ggsave(file.path(output_dir, "PVL_recovery_lambda_sd.png"), lambda_sd_plots, 
       width = 12, height = 10, dpi = 300)
cat("Saved: PVL_recovery_lambda_sd.png\n")


# Save results to RData file
save(true_mu_w, true_mu_A, true_mu_theta, true_mu_a,
     infer_mu_w, infer_mu_A, infer_mu_theta, infer_mu_a,
     true_lambda_w, true_lambda_A, true_lambda_theta, true_lambda_a,
     infer_lambda_w, infer_lambda_A, infer_lambda_theta, infer_lambda_a,
     file = file.path(output_dir, "PVL_recovery_results.RData"))
cat("Saved: PVL_recovery_results.RData\n")


# Calculate recovery metrics and save to text file
recovery_file <- file.path(output_dir, "recovery_metrics.txt")
sink(recovery_file)

cat("\n=== RECOVERY METRICS (50 iterations) ===\n")

cat("\nCorrelations (need > 0.7 for good recovery):\n")
cat("w:    ", round(cor(true_mu_w, infer_mu_w), 3), "\n")
cat("A:    ", round(cor(true_mu_A, infer_mu_A), 3), "\n")
cat("theta:", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("a:    ", round(cor(true_mu_a, infer_mu_a), 3), "\n")

cat("\nBias (should be near 0):\n")
cat("w:    ", round(mean(infer_mu_w - true_mu_w), 3), "\n")
cat("A:    ", round(mean(infer_mu_A - true_mu_A), 3), "\n")
cat("theta:", round(mean(infer_mu_theta - true_mu_theta), 3), "\n")
cat("a:    ", round(mean(infer_mu_a - true_mu_a), 3), "\n")

cat("\nRMSE (lower is better):\n")
cat("w:    ", round(sqrt(mean((infer_mu_w - true_mu_w)^2)), 3), "\n")
cat("A:    ", round(sqrt(mean((infer_mu_A - true_mu_A)^2)), 3), "\n")
cat("theta:", round(sqrt(mean((infer_mu_theta - true_mu_theta)^2)), 3), "\n")
cat("a:    ", round(sqrt(mean((infer_mu_a - true_mu_a)^2)), 3), "\n")

sink()
cat("Saved: recovery_metrics.txt\n")
