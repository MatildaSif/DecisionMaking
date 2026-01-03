install.packages("pacman")
install.packages("truncnorm")
pacman::p_load(extraDistr, R2jags, parallel, ggpubr, truncnorm)

set.seed(1702)

setwd("/work/JoMat/DecisionMaking/src/VSE")

# defining a function for calculating the maximum of the posterior density
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

# Simple recovery plot function
recov_plot <- function(true_vals, infer_vals, labels, title) {
  df <- data.frame(true = true_vals, infer = infer_vals)
  
  ggplot(df, aes(x = true, y = infer)) +
    geom_point(size = 2, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(x = labels[1], y = labels[2], title = title) +
    theme_minimal() +
    theme(aspect.ratio = 1)
}

#------ create task environment -------------------
# Bad frequent
A_R <- rep(100,10)
A_L <- c(rep(-250,5),rep(0,5))
# Bad infrequent
B_R <- rep(100,10)
B_L <- c(rep(-1250,1),rep(0,9))
# Good frequent
C_R <- rep(50,10)
C_L <- c(rep(-50,5),rep(0,5))
# Good infrequent
D_R <- rep(50,10)
D_L <- c(rep(-250,1),rep(0,9))

A <- c()
for (i in 1:10) {A <- (append(A,A_R + sample(A_L)))}
B <- c()
for (i in 1:10) {B <- (append(B,B_R + sample(B_L)))}
C <- c()
for (i in 1:10) {C <- (append(C,C_R + sample(C_L)))}
D <- c()
for (i in 1:10) {D <- (append(D,D_R + sample(D_L)))}

# extract rewards from reward trials
RA <- A
RA[A < 0] = 0
RB <- B
RB[B < 0] = 0
RC <- C
RC[C < 0] = 0
RD <- D
RD[D < 0] = 0

# extract losses from reward trials
LA <- A
LA[A > 0] <- 0
LA <- abs(LA)
LB <- B
LB[B > 0] <- 0
LB <- abs(LB)
LC <- C
LC[C > 0] <- 0
LC <- abs(LC)
LD <- D
LD[D > 0] <- 0
LD <- abs(LD)

# Create payoff matrices (deck-specific)
R_decks <- cbind(RA, RB, RC, RD) # We have 100 rows and 4 decs (columns)
L_decks <- cbind(LA, LB, LC, LD) # We have the same here

###--------------Run full parameter recovery -------------
ntrials <- 100 
niterations <- 2 # number of recovery iterations (use 50 for full recovery)
nsubs <- 48 # number of simulated subjects
ntrials_all <- rep(ntrials, nsubs) # all 48 simulated subs have 100 trials each

# Initialize arrays for GROUP-LEVEL parameters (mu = mean, lambda = precision)
# mu (group means)
true_mu_theta <- array(NA, c(niterations))
true_mu_delta <- array(NA, c(niterations))
true_mu_alpha <- array(NA, c(niterations))
true_mu_phi <- array(NA, c(niterations))
true_mu_c <- array(NA, c(niterations))

infer_mu_theta <- array(NA, c(niterations))
infer_mu_delta <- array(NA, c(niterations))
infer_mu_alpha <- array(NA, c(niterations))
infer_mu_phi <- array(NA, c(niterations))
infer_mu_c <- array(NA, c(niterations))

# sigma (between-subject SD) - storing the TRUE sigma values
true_sigma_theta <- array(NA, c(niterations))
true_sigma_delta <- array(NA, c(niterations))
true_sigma_alpha <- array(NA, c(niterations))
true_sigma_phi <- array(NA, c(niterations))
true_sigma_c <- array(NA, c(niterations))

# lambda (precision for JAGS) - inferred values
infer_lambda_theta <- array(NA, c(niterations))
infer_lambda_delta <- array(NA, c(niterations))
infer_lambda_alpha <- array(NA, c(niterations))
infer_lambda_phi <- array(NA, c(niterations))
infer_lambda_c <- array(NA, c(niterations))

# Source the simulation function ONCE before the loop
source('heir_VSE_sim.R')

start_time = Sys.time()

for (i in 1:niterations) {
  
  # Sample GROUP-LEVEL parameters (hyperparameters)
  mu_theta <- runif(1, 0.1, 0.5)      # learning rate (0-1)
  mu_delta <- runif(1, 0.5, 2.0)      # exploration bonus (positive)
  mu_alpha <- runif(1, 0.1, 0.5)      # exploration weight (0-1)
  mu_phi <- runif(1, 0.1, 0.5)        # decay rate (0-1)
  mu_c <- runif(1, 0.5, 2.0)          # inverse temperature (positive)
  
  # Between-subject variability (sigma)
  sigma_theta <- runif(1, 0.05, 0.15)
  sigma_delta <- runif(1, 0.1, 0.4)
  sigma_alpha <- runif(1, 0.05, 0.15)
  sigma_phi <- runif(1, 0.05, 0.15)
  sigma_c <- runif(1, 0.1, 0.3)
  
  # Run simulation (R_decks and L_decks are deck-specific payoff matrices)
  VSE_sims <- heir_VSE_sim(R_decks, L_decks, nsubs, ntrials_all,
                           mu_theta, mu_delta, mu_alpha, mu_phi, mu_c,
                           sigma_theta, sigma_delta, sigma_alpha, sigma_phi, sigma_c)
  
  
  x <- VSE_sims$x
  R <- VSE_sims$R  # realized rewards (nsubs x ntrials)
  L <- VSE_sims$L  # realized losses (nsubs x ntrials)
  
  print(dim(R))
  print(dim(L))
  
  # Set up JAGS
  data <- list("x", "R", "L", "ntrials_all", "nsubs")
  params <- c("mu_theta", "mu_delta", "mu_alpha", "mu_phi", "mu_c",
              "lambda_theta", "lambda_delta", "lambda_alpha", "lambda_phi", "lambda_c")
  
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file="heir_VSE.txt", n.chains=3, 
                           n.iter=1000, n.burnin=500, n.thin=1, n.cluster=4)
  
  # Store true GROUP-LEVEL parameters
  true_mu_theta[i] <- mu_theta
  true_mu_delta[i] <- mu_delta
  true_mu_alpha[i] <- mu_alpha
  true_mu_phi[i] <- mu_phi
  true_mu_c[i] <- mu_c
  
  # Store true sigma values
  true_sigma_theta[i] <- sigma_theta
  true_sigma_delta[i] <- sigma_delta
  true_sigma_alpha[i] <- sigma_alpha
  true_sigma_phi[i] <- sigma_phi
  true_sigma_c[i] <- sigma_c
  
  # Extract inferred GROUP-LEVEL parameters (MPD)
  Y <- samples$BUGSoutput$sims.list
  
  infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_delta[i] <- MPD(Y$mu_delta)
  infer_mu_alpha[i] <- MPD(Y$mu_alpha)
  infer_mu_phi[i] <- MPD(Y$mu_phi)
  infer_mu_c[i] <- MPD(Y$mu_c)
  
  infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_delta[i] <- MPD(Y$lambda_delta)
  infer_lambda_alpha[i] <- MPD(Y$lambda_alpha)
  infer_lambda_phi[i] <- MPD(Y$lambda_phi)
  infer_lambda_c[i] <- MPD(Y$lambda_c)
  
  print(paste("Iteration", i, "completed"))
}

end_time = Sys.time()
print(end_time - start_time)

# Convert lambda to sigma for comparison (sigma = 1/sqrt(lambda))
infer_sigma_theta <- 1/sqrt(infer_lambda_theta)
infer_sigma_delta <- 1/sqrt(infer_lambda_delta)
infer_sigma_alpha <- 1/sqrt(infer_lambda_alpha)
infer_sigma_phi <- 1/sqrt(infer_lambda_phi)
infer_sigma_c <- 1/sqrt(infer_lambda_c)

# Plotting recovery results for mu parameters
pl1 <- recov_plot(true_mu_theta, infer_mu_theta, 
                  c("true mu_theta", "infer mu_theta"), 'mu_theta recovery')
pl2 <- recov_plot(true_mu_delta, infer_mu_delta, 
                  c("true mu_delta", "infer mu_delta"), 'mu_delta recovery')
pl3 <- recov_plot(true_mu_alpha, infer_mu_alpha, 
                  c("true mu_alpha", "infer mu_alpha"), 'mu_alpha recovery')
pl4 <- recov_plot(true_mu_phi, infer_mu_phi, 
                  c("true mu_phi", "infer mu_phi"), 'mu_phi recovery')
pl5 <- recov_plot(true_mu_c, infer_mu_c, 
                  c("true mu_c", "infer mu_c"), 'mu_c recovery')

# Arrange and display
mu_plot <- ggarrange(pl1, pl2, pl3, pl4, pl5, ncol = 3, nrow = 2)
print(mu_plot)

# Save mu parameters plot
ggsave("VSE_recovery_mu_parameters.png", mu_plot, 
       width = 12, height = 8, dpi = 300)
cat("\nSaved: VSE_recovery_mu_parameters.png\n")


# Plotting recovery results for sigma parameters
pl6 <- recov_plot(true_sigma_theta, infer_sigma_theta, 
                  c("true sigma_theta", "infer sigma_theta"), 'sigma_theta recovery')
pl7 <- recov_plot(true_sigma_delta, infer_sigma_delta, 
                  c("true sigma_delta", "infer sigma_delta"), 'sigma_delta recovery')
pl8 <- recov_plot(true_sigma_alpha, infer_sigma_alpha, 
                  c("true sigma_alpha", "infer sigma_alpha"), 'sigma_alpha recovery')
pl9 <- recov_plot(true_sigma_phi, infer_sigma_phi, 
                  c("true sigma_phi", "infer sigma_phi"), 'sigma_phi recovery')
pl10 <- recov_plot(true_sigma_c, infer_sigma_c, 
                   c("true sigma_c", "infer sigma_c"), 'sigma_c recovery')

# Arrange and display
sigma_plot <- ggarrange(pl6, pl7, pl8, pl9, pl10, ncol = 3, nrow = 2)
print(sigma_plot)

# Save sigma parameters plot
ggsave("VSE_recovery_sigma_parameters.png", sigma_plot, 
       width = 12, height = 8, dpi = 300)
cat("Saved: VSE_recovery_sigma_parameters.png\n")

# Save results
save(true_mu_theta, true_mu_delta, true_mu_alpha, true_mu_phi, true_mu_c,
     infer_mu_theta, infer_mu_delta, infer_mu_alpha, infer_mu_phi, infer_mu_c,
     true_sigma_theta, true_sigma_delta, true_sigma_alpha, true_sigma_phi, true_sigma_c,
     infer_sigma_theta, infer_sigma_delta, infer_sigma_alpha, infer_sigma_phi, infer_sigma_c,
     infer_lambda_theta, infer_lambda_delta, infer_lambda_alpha, infer_lambda_phi, infer_lambda_c,
     file = "VSE_recovery_results.RData")

# Print summary statistics
cat("\n=== RECOVERY SUMMARY ===\n")
cat("\nGroup Means (mu) - Correlations:\n")
cat("theta:", cor(true_mu_theta, infer_mu_theta), "\n")
cat("delta:", cor(true_mu_delta, infer_mu_delta), "\n")
cat("alpha:", cor(true_mu_alpha, infer_mu_alpha), "\n")
cat("phi:", cor(true_mu_phi, infer_mu_phi), "\n")
cat("c:", cor(true_mu_c, infer_mu_c), "\n")

cat("\nGroup SDs (sigma) - Correlations:\n")
cat("theta:", cor(true_sigma_theta, infer_sigma_theta), "\n")
cat("delta:", cor(true_sigma_delta, infer_sigma_delta), "\n")
cat("alpha:", cor(true_sigma_alpha, infer_sigma_alpha), "\n")
cat("phi:", cor(true_sigma_phi, infer_sigma_phi), "\n")
cat("c:", cor(true_sigma_c, infer_sigma_c), "\n")