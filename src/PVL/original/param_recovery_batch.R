#!/usr/bin/env Rscript

# Get seed/iteration from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  start_iter <- 1
  seed_offset <- 0
} else {
  seed_offset <- as.numeric(args[1])
  start_iter <- 1
}

# Set seed based on offset
set.seed(1983 + seed_offset)

pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)

setwd('/work/JoMat/DecisionMaking/src/PVL/original/')

# defining a function for calculating the maximum of the posterior density
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
ntrials <- 100
nstruct <- 10
freq <- 0.5
infreq <- 0.1
bad_r <- 100
bad_freq_l <- -250
bad_infreq_l <- -1250
good_r <- 50
good_freq_l <- -50
good_infreq_l <- -250

# Bad frequent
A_R <- rep(bad_r, nstruct)
A_L <- c(rep(bad_freq_l, nstruct*freq), rep(0, nstruct*(1-freq)))

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq)))

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq), rep(0, nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq)))

# Create pseudorandomized payoff structure
A <- array(NA, ntrials)
B <- array(NA, ntrials)
C <- array(NA, ntrials)
D <- array(NA, ntrials)

for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L))
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}

payoff <- cbind(A, B, C, D) / 100

###----- Run parameter recovery for this iteration batch -----
# Running 10 iterations per job for reasonable runtime
niterations <- 10
nsubs <- 48
ntrials_all <- rep(100, 48)

# Initialize storage arrays
true_mu_w <- array(NA, c(niterations))
true_mu_A <- array(NA, c(niterations))
true_mu_theta <- array(NA, c(niterations))
true_mu_a <- array(NA, c(niterations))

infer_mu_w <- array(NA, c(niterations))
infer_mu_A <- array(NA, c(niterations))
infer_mu_theta <- array(NA, c(niterations))
infer_mu_a <- array(NA, c(niterations))

true_lambda_w <- array(NA, c(niterations))
true_lambda_A <- array(NA, c(niterations))
true_lambda_theta <- array(NA, c(niterations))
true_lambda_a <- array(NA, c(niterations))

infer_lambda_w <- array(NA, c(niterations))
infer_lambda_A <- array(NA, c(niterations))
infer_lambda_theta <- array(NA, c(niterations))
infer_lambda_a <- array(NA, c(niterations))

cat("Starting parameter recovery batch with seed offset:", seed_offset, "\n")
start_time <- Sys.time()

for (i in 1:niterations) {
  ntrials <- ntrials_all
  
  mu_w <- runif(1, .5, 2.5)
  mu_A <- runif(1, 0, 1)
  mu_theta <- runif(1, 0, 2)
  mu_a <- runif(1, 0, 1)
  
  sigma_w <- runif(1, 0, 0.2)
  sigma_A <- runif(1, 0, 0.1)
  sigma_theta <- runif(1, 0, 0.2)
  sigma_a <- runif(1, 0, 0.1)
  
  source('hier_PVL_sim.R')
  PVL_sims <- hier_PVL_sim(payoff, nsubs, ntrials, mu_w, mu_A, mu_a, mu_theta,
                           sigma_w, sigma_A, sigma_a, sigma_theta)
  
  x <- PVL_sims$x
  X <- PVL_sims$X
  
  # Run JAGS model
  data <- list("x", "X", "ntrials", "nsubs")
  params <- c("mu_w", "mu_A", "mu_theta", "mu_a", "lambda_w", "lambda_A", "lambda_theta", "lambda_a")
  samples <- jags.parallel(data, inits = NULL, params,
                           model.file = "hier_PVL.txt", n.chains = 3,
                           n.iter = 3000, n.burnin = 1000, n.thin = 1, n.cluster = 4)
  
  # Store true values
  true_mu_w[i] <- mu_w
  true_mu_A[i] <- mu_A
  true_mu_theta[i] <- mu_theta
  true_mu_a[i] <- mu_a
  
  # Extract MAP estimates
  Y <- samples$BUGSoutput$sims.list
  infer_mu_w[i] <- MPD(Y$mu_w)
  infer_mu_A[i] <- MPD(Y$mu_A)
  infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_a[i] <- MPD(Y$mu_a)
  
  # Store true lambda/sigma values
  true_lambda_w[i] <- sigma_w
  true_lambda_A[i] <- sigma_A
  true_lambda_theta[i] <- sigma_theta
  true_lambda_a[i] <- sigma_a
  
  # Extract MAP estimates for precision
  infer_lambda_w[i] <- MPD(Y$lambda_w)
  infer_lambda_A[i] <- MPD(Y$lambda_A)
  infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_a[i] <- MPD(Y$lambda_a)
  
  cat(sprintf("Batch %d, Iteration %d complete\n", seed_offset, i))
}

end_time <- Sys.time()
cat("Time elapsed:", format(end_time - start_time), "\n")

# Save results to disk
output_file <- sprintf('results/param_recovery_batch_%03d.RData', seed_offset)
dir.create('results', showWarnings = FALSE)
save(true_mu_w, true_mu_A, true_mu_theta, true_mu_a,
     infer_mu_w, infer_mu_A, infer_mu_theta, infer_mu_a,
     true_lambda_w, true_lambda_A, true_lambda_theta, true_lambda_a,
     infer_lambda_w, infer_lambda_A, infer_lambda_theta, infer_lambda_a,
     file = output_file)

cat("Results saved to:", output_file, "\n")