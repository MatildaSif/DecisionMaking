install.packages("pacman")
pacman::p_load(extraDistr, R2jags, parallel, ggpubr, rstatix, gridExtra, polspline, truncnorm, ggplot2, tidyr, dplyr)

set.seed(1702)

setwd('/work/JoMat/DecisionMaking/src/VSE')

# defining a function for calculating the maximum of the posterior density
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------load data----
ctr_data <- read.table("../../data/IGTdata_healthy_control_clean.txt", header=TRUE)
opi_data <- read.table("../../data/IGTdata_heroin_clean.txt", header=TRUE)
amp_data <- read.table("../../data/IGTdata_amphetamine_clean.txt", header=TRUE)

#----------prepare data for jags models----
# identify and count unique subject IDs
subIDs_ctr <- unique(ctr_data$subjID)
nsubs_ctr <- length(subIDs_ctr)

subIDs_opi <- unique(opi_data$subjID)
nsubs_opi <- length(subIDs_opi)

subIDs_amp <- unique(amp_data$subjID)
nsubs_amp <- length(subIDs_amp)

ntrials_max <- 100

# Get raw data
x_raw_ctr <- ctr_data$deck
X_raw_ctr <- (ctr_data$gain + ctr_data$loss)

x_raw_opi <- opi_data$deck
X_raw_opi <- (opi_data$gain + opi_data$loss)

x_raw_amp <- amp_data$deck
X_raw_amp <- (amp_data$gain + amp_data$loss)

#--- Assign choices and outcomes in trial x sub matrix ---

# Empty arrays to fill
ntrials_ctr <- array(0, c(nsubs_ctr))
x_ctr <- array(0, c(nsubs_ctr, ntrials_max))
X_ctr <- array(0, c(nsubs_ctr, ntrials_max))

ntrials_opi <- array(0, c(nsubs_opi))
x_opi <- array(0, c(nsubs_opi, ntrials_max))
X_opi <- array(0, c(nsubs_opi, ntrials_max))

ntrials_amp <- array(0, c(nsubs_amp))
x_amp <- array(0, c(nsubs_amp, ntrials_max))
X_amp <- array(0, c(nsubs_amp, ntrials_max))

# Make control data matrices
for (s in 1:nsubs_ctr) {
  
  # Record n trials for subject s
  ntrials_ctr[s] <- length(x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]])
  
  # Pad trials with NA if n trials < maximum
  x_sub_ctr <- x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(x_sub_ctr) <- ntrials_max
  
  X_sub_ctr <- X_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(X_sub_ctr) <- ntrials_max
  
  # Assign arrays
  x_ctr[s,] <- x_sub_ctr
  X_ctr[s,] <- X_sub_ctr
}

# Make heroin data matrices
for (s in 1:nsubs_opi) {
  
  ntrials_opi[s] <- length(x_raw_opi[opi_data$subjID == subIDs_opi[s]])
  
  x_sub_opi <- x_raw_opi[opi_data$subjID == subIDs_opi[s]] 
  length(x_sub_opi) <- ntrials_max
  
  X_sub_opi <- X_raw_opi[opi_data$subjID == subIDs_opi[s]] 
  length(X_sub_opi) <- ntrials_max
  
  x_opi[s,] <- x_sub_opi
  X_opi[s,] <- X_sub_opi
}

# Make amphetamine data matrices
for (s in 1:nsubs_amp) {
  
  ntrials_amp[s] <- length(x_raw_amp[amp_data$subjID == subIDs_amp[s]])
  
  x_sub_amp <- x_raw_amp[amp_data$subjID == subIDs_amp[s]] 
  length(x_sub_amp) <- ntrials_max
  
  X_sub_amp <- X_raw_amp[amp_data$subjID == subIDs_amp[s]] 
  length(X_sub_amp) <- ntrials_max
  
  x_amp[s,] <- x_sub_amp
  X_amp[s,] <- X_sub_amp
}

# ***  TEST MODE - so I can built it hehe ***
TEST_MODE <- F  # Set to FALSE for full run

if (TEST_MODE) {
  cat("\n*** RUNNING IN TEST MODE ***\n")
  
  # Store original counts
  nsubs_ctr_full <- nsubs_ctr
  nsubs_opi_full <- nsubs_opi
  nsubs_amp_full <- nsubs_amp
  
  # Reduce sample size
  nsubs_ctr <- min(3, nsubs_ctr)
  nsubs_opi <- min(3, nsubs_opi)
  nsubs_amp <- min(3, nsubs_amp)
  
  cat("Reducing from", nsubs_ctr_full, "to", nsubs_ctr, "controls\n")
  cat("Reducing from", nsubs_opi_full, "to", nsubs_opi, "heroin\n")
  cat("Reducing from", nsubs_amp_full, "to", nsubs_amp, "amphetamine\n")
  
  # Subset the data
  x_ctr <- x_ctr[1:nsubs_ctr, ]
  X_ctr <- X_ctr[1:nsubs_ctr, ]
  ntrials_ctr <- ntrials_ctr[1:nsubs_ctr]
  
  x_opi <- x_opi[1:nsubs_opi, ]
  X_opi <- X_opi[1:nsubs_opi, ]
  ntrials_opi <- ntrials_opi[1:nsubs_opi]
  
  x_amp <- x_amp[1:nsubs_amp, ]
  X_amp <- X_amp[1:nsubs_amp, ]
  ntrials_amp <- ntrials_amp[1:nsubs_amp]
  
}

#----------Run VSE model comparison----
cat("\n=== Running VSE Group Comparison Model ===\n")
cat("Controls:", nsubs_ctr, "subjects\n")
cat("Heroin:", nsubs_opi, "subjects\n")
cat("Amphetamine:", nsubs_amp, "subjects\n\n")

# Set up JAGS data
data <- list("x_ctr", "X_ctr", "ntrials_ctr", "nsubs_ctr",
             "x_opi", "X_opi", "ntrials_opi", "nsubs_opi",
             "x_amp", "X_amp", "ntrials_amp", "nsubs_amp")

# Parameters to track - INCLUDING p for posterior predictive checks
params <- c("mu_theta_ctr", "mu_delta_ctr", "mu_alpha_ctr", "mu_phi_ctr", "mu_c_ctr",
            "mu_theta_opi", "mu_delta_opi", "mu_alpha_opi", "mu_phi_opi", "mu_c_opi",
            "mu_theta_amp", "mu_delta_amp", "mu_alpha_amp", "mu_phi_amp", "mu_c_amp",
            "diff_theta_ctr_opi", "diff_delta_ctr_opi", "diff_alpha_ctr_opi", 
            "diff_phi_ctr_opi", "diff_c_ctr_opi",
            "diff_theta_ctr_amp", "diff_delta_ctr_amp", "diff_alpha_ctr_amp",
            "diff_phi_ctr_amp", "diff_c_ctr_amp",
            "diff_theta_opi_amp", "diff_delta_opi_amp", "diff_alpha_opi_amp",
            "diff_phi_opi_amp", "diff_c_opi_amp",
            "theta_ctr", "delta_ctr", "alpha_ctr", "phi_ctr", "c_ctr",
            "theta_opi", "delta_opi", "alpha_opi", "phi_opi", "c_opi",
            "theta_amp", "delta_amp", "alpha_amp", "phi_amp", "c_amp")  


# Run JAGS
start_time <- Sys.time()

if (TEST_MODE) {
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file="VSE_group_compare.txt",
                           n.chains=2, 
                           n.iter=100, 
                           n.burnin=50, 
                           n.thin=1,
                           n.cluster=2)
} else {
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file="VSE_group_compare.txt",
                           n.chains=3, 
                           n.iter=5000, 
                           n.burnin=1500, 
                           n.thin=1,
                           n.cluster=12)
}
end_time <- Sys.time()

cat("\nModel fitting completed in:", 
    round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")

#----------Extract and summarize results----
Y <- samples$BUGSoutput$sims.list

# Group-level parameters
cat("\n=== GROUP-LEVEL PARAMETER ESTIMATES ===\n")

cat("\n--- CONTROLS ---\n")
cat("theta:", round(MPD(Y$mu_theta_ctr), 3), "\n")
cat("delta:", round(MPD(Y$mu_delta_ctr), 3), "\n")
cat("alpha:", round(MPD(Y$mu_alpha_ctr), 3), "\n")
cat("phi:  ", round(MPD(Y$mu_phi_ctr), 3), "\n")
cat("c:    ", round(MPD(Y$mu_c_ctr), 3), "\n")

cat("\n--- HEROIN USERS ---\n")
cat("theta:", round(MPD(Y$mu_theta_opi), 3), "\n")
cat("delta:", round(MPD(Y$mu_delta_opi), 3), "\n")
cat("alpha:", round(MPD(Y$mu_alpha_opi), 3), "\n")
cat("phi:  ", round(MPD(Y$mu_phi_opi), 3), "\n")
cat("c:    ", round(MPD(Y$mu_c_opi), 3), "\n")

cat("\n--- AMPHETAMINE USERS ---\n")
cat("theta:", round(MPD(Y$mu_theta_amp), 3), "\n")
cat("delta:", round(MPD(Y$mu_delta_amp), 3), "\n")
cat("alpha:", round(MPD(Y$mu_alpha_amp), 3), "\n")
cat("phi:  ", round(MPD(Y$mu_phi_amp), 3), "\n")
cat("c:    ", round(MPD(Y$mu_c_amp), 3), "\n")

#----------Visualize group differences----

# Function to create comparison plot
compare_plot <- function(param_name, ctr, opi, amp, xlim_range) {
  plot(density(ctr), xlim=xlim_range, main=param_name, 
       col="lightgreen", lwd=2, xlab="Parameter value", ylab="Density")
  lines(density(opi), col="lightblue", lwd=2)
  lines(density(amp), col="salmon", lwd=2)
  legend("topright", legend=c("Controls", "Heroin", "Amphetamine"), 
         col=c("lightgreen", "lightblue", "salmon"), lwd=2)
  abline(v=0, lty=2, col="black")
}


# Posterior plots 
#----------Prior vs Posterior Plots----

# Function to plot prior vs posterior
plot_prior_posterior <- function(posterior_samples, prior_mean, prior_sd, 
                                 lower_bound = -Inf, upper_bound = Inf,
                                 param_name, xlim_range) {
  
  # Create prior distribution
  x_seq <- seq(xlim_range[1], xlim_range[2], length.out = 1000)
  
  if (is.finite(lower_bound) && is.finite(upper_bound)) {
    # Truncated normal
    prior_density <- dtruncnorm(x_seq, a = lower_bound, b = upper_bound, 
                                mean = prior_mean, sd = prior_sd)
  } else if (is.finite(lower_bound)) {
    # Lower truncated
    prior_density <- dtruncnorm(x_seq, a = lower_bound, b = Inf, 
                                mean = prior_mean, sd = prior_sd)
  } else {
    # Normal (unbounded)
    prior_density <- dnorm(x_seq, mean = prior_mean, sd = prior_sd)
  }
  
  # Plot posterior (filled)
  plot(density(posterior_samples), main = param_name, 
       xlab = "Parameter value", ylab = "Density",
       col = "red", lwd = 2, xlim = xlim_range, ylim = c(0, max(density(posterior_samples)$y) * 1.1))
  polygon(density(posterior_samples), col = rgb(0, 0, 1, 0.3), border = NA)
  
  # Add prior (dashed line)
  lines(x_seq, prior_density, col = "black", lwd = 2, lty = 2)
  
  legend("topright", legend = c("Posterior", "Prior"), 
         col = c("red", "black"), lwd = 2, lty = c(1, 2), bty = "n")
}

# Create prior-posterior plots for CONTROLS
png("VSE_prior_posterior_controls.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

plot_prior_posterior(Y$mu_theta_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = 1,
                     param_name = "theta (Controls)", xlim_range = c(0, 1))

plot_prior_posterior(Y$mu_delta_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = 1,
                     param_name = "delta (Controls)", xlim_range = c(0, 1))

plot_prior_posterior(Y$mu_alpha_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = 1,
                     param_name = "alpha (Controls)", xlim_range = c(0, 1))

plot_prior_posterior(Y$mu_phi_ctr, prior_mean = 0, prior_sd = sqrt(1/0.1), 
                     lower_bound = -Inf, upper_bound = Inf,
                     param_name = "phi (Controls)", xlim_range = c(-3, 3))

plot_prior_posterior(Y$mu_c_ctr, prior_mean = 0, prior_sd = sqrt(1/0.1), 
                     lower_bound = 0, upper_bound = Inf,
                     param_name = "c (Controls)", xlim_range = c(0, 10))

dev.off()
cat("Saved: VSE_prior_posterior_controls.png\n")

# Repeat for HEROIN and AMPHETAMINE groups
png("VSE_prior_posterior_heroin.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

plot_prior_posterior(Y$mu_theta_opi, 0, 1, 0, 1, "theta (Heroin)", c(0, 1))
plot_prior_posterior(Y$mu_delta_opi, 0, 1, 0, 1, "delta (Heroin)", c(0, 1))
plot_prior_posterior(Y$mu_alpha_opi, 0, 1, 0, 1, "alpha (Heroin)", c(0, 1))
plot_prior_posterior(Y$mu_phi_opi, 0, sqrt(1/0.1), -Inf, Inf, "phi (Heroin)", c(-3, 3))
plot_prior_posterior(Y$mu_c_opi, 0, sqrt(1/0.1), 0, Inf, "c (Heroin)", c(0, 10))

dev.off()
cat("Saved: VSE_prior_posterior_heroin.png\n")

png("VSE_prior_posterior_amphetamine.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

plot_prior_posterior(Y$mu_theta_amp, 0, 1, 0, 1, "theta (Amphetamine)", c(0, 1))
plot_prior_posterior(Y$mu_delta_amp, 0, 1, 0, 1, "delta (Amphetamine)", c(0, 1))
plot_prior_posterior(Y$mu_alpha_amp, 0, 1, 0, 1, "alpha (Amphetamine)", c(0, 1))
plot_prior_posterior(Y$mu_phi_amp, 0, sqrt(1/0.1), -Inf, Inf, "phi (Amphetamine)", c(-3, 3))
plot_prior_posterior(Y$mu_c_amp, 0, sqrt(1/0.1), 0, Inf, "c (Amphetamine)", c(0, 10))

dev.off()
cat("Saved: VSE_prior_posterior_amphetamine.png\n")

# ============================================================
# Sequential Exploration (SE3 / SE4) – ggplot2 version
# ============================================================


# ---- Helper functions ----

all_unique <- function(x) {
  x <- na.omit(x)
  length(unique(x)) == length(x)
}

compute_SE <- function(choices, ntrials, ntrials_max) {
  
  SE3 <- rep(NA, ntrials_max)
  SE4 <- rep(NA, ntrials_max)
  
  # SE3: rolling window of length 3 → assign to trial t
  if (ntrials >= 3) {
    for (t in 3:ntrials) {
      window <- choices[(t-2):t]
      if (sum(!is.na(window)) == 3) {
        SE3[t] <- as.numeric(all_unique(window))
      }
    }
  }
  
  # SE4: rolling window of length 4 → assign to trial t
  if (ntrials >= 4) {
    for (t in 4:ntrials) {
      window <- choices[(t-3):t]
      if (sum(!is.na(window)) == 4) {
        SE4[t] <- as.numeric(all_unique(window))
      }
    }
  }
  
  list(SE3 = SE3, SE4 = SE4)
}

compute_group_SE <- function(x_mat, ntrials_vec, nsubs, ntrials_max) {
  
  SE3_mat <- matrix(NA, nsubs, ntrials_max)
  SE4_mat <- matrix(NA, nsubs, ntrials_max)
  
  for (s in 1:nsubs) {
    se <- compute_SE(x_mat[s, ], ntrials_vec[s], ntrials_max)  # Added ntrials_max
    SE3_mat[s, ] <- se$SE3
    SE4_mat[s, ] <- se$SE4
  }
  
  list(SE3 = SE3_mat, SE4 = SE4_mat)
}

# Calculate mean across subjects for each trial
summarize_SE <- function(SE_mat) {
  
  mean_SE <- colMeans(SE_mat, na.rm = TRUE)
  n_eff <- colSums(!is.na(SE_mat))
  se_SE <- apply(SE_mat, 2, sd, na.rm = TRUE) / sqrt(n_eff)
  
  list(
    mean  = mean_SE,
    lower = mean_SE - 1.96 * se_SE,
    upper = mean_SE + 1.96 * se_SE
  )
}


# ---- Compute SE per group ----

SE_ctr <- compute_group_SE(x_ctr, ntrials_ctr, nsubs_ctr, ntrials_max)
SE_opi <- compute_group_SE(x_opi, ntrials_opi, nsubs_opi, ntrials_max)
SE_amp <- compute_group_SE(x_amp, ntrials_amp, nsubs_amp, ntrials_max)

SE3_ctr_sum <- summarize_SE(SE_ctr$SE3)
SE3_opi_sum <- summarize_SE(SE_opi$SE3)
SE3_amp_sum <- summarize_SE(SE_amp$SE3)

SE4_ctr_sum <- summarize_SE(SE_ctr$SE4)
SE4_opi_sum <- summarize_SE(SE_opi$SE4)
SE4_amp_sum <- summarize_SE(SE_amp$SE4)

# ---- Prepare data for ggplot ----

# Create dataframe for SE3
df_SE3 <- data.frame(
  trial = rep(1:ntrials_max, 3),
  mean = c(SE3_ctr_sum$mean, SE3_opi_sum$mean, SE3_amp_sum$mean),
  lower = c(SE3_ctr_sum$lower, SE3_opi_sum$lower, SE3_amp_sum$lower),
  upper = c(SE3_ctr_sum$upper, SE3_opi_sum$upper, SE3_amp_sum$upper),
  group = rep(c("Controls", "Heroin", "Amphetamine"), each = ntrials_max)
)

# Create dataframe for SE4
df_SE4 <- data.frame(
  trial = rep(1:ntrials_max, 3),
  mean = c(SE4_ctr_sum$mean, SE4_opi_sum$mean, SE4_amp_sum$mean),
  lower = c(SE4_ctr_sum$lower, SE4_opi_sum$lower, SE4_amp_sum$lower),
  upper = c(SE4_ctr_sum$upper, SE4_opi_sum$upper, SE4_amp_sum$upper),
  group = rep(c("Controls", "Heroin", "Amphetamine"), each = ntrials_max)
)


# Create dataframe for SE3
df_SE3 <- data.frame(
  trial = rep(1:ntrials_max, 3),
  mean = c(SE3_ctr_sum$mean, SE3_opi_sum$mean, SE3_amp_sum$mean),
  lower = c(SE3_ctr_sum$lower, SE3_opi_sum$lower, SE3_amp_sum$lower),
  upper = c(SE3_ctr_sum$upper, SE3_opi_sum$upper, SE3_amp_sum$upper),
  group = rep(c("Controls", "Heroin", "Amphetamine"), each = ntrials_max)
)

# Create dataframe for SE4
df_SE4 <- data.frame(
  trial = rep(1:ntrials_max, 3),
  mean = c(SE4_ctr_sum$mean, SE4_opi_sum$mean, SE4_amp_sum$mean),
  lower = c(SE4_ctr_sum$lower, SE4_opi_sum$lower, SE4_amp_sum$lower),
  upper = c(SE4_ctr_sum$upper, SE4_opi_sum$upper, SE4_amp_sum$upper),
  group = rep(c("Controls", "Heroin", "Amphetamine"), each = ntrials_max)
)

# ---- Plot SE3 ----

p_SE3 <- ggplot(df_SE3, aes(x = trial, y = mean, color = group, fill = group)) +
  geom_smooth(method = "loess", span = 0.15, se = TRUE, alpha = 0.2, size = 1.2) +
  geom_hline(yintercept = 0.33, linetype = "dashed", color = "gray30", size = 0.8) +
  scale_color_manual(values = c("Controls" = "lightgreen", "Heroin" = "lightblue", "Amphetamine" = "salmon")) +
  scale_fill_manual(values = c("Controls" = "lightgreen", "Heroin" = "lightblue", "Amphetamine" = "salmon")) +
  labs(
    title = "SE3 across IGT",
    x = "Trial",
    y = "Frequency",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 0.85), xlim = c(0, 100))

# ---- Plot SE4 ----

p_SE4 <- ggplot(df_SE4, aes(x = trial, y = mean, color = group, fill = group)) +
  geom_smooth(method = "loess", span = 0.15, se = TRUE, alpha = 0.2, size = 1.2) +
  geom_hline(yintercept = 0.09, linetype = "dashed", color = "gray30", size = 0.8) +
  scale_color_manual(values = c("Controls" = "lightgreen", "Heroin" = "lightblue", "Amphetamine" = "salmon")) +
  scale_fill_manual(values = c("Controls" = "lightgreen", "Heroin" = "lightblue", "Amphetamine" = "salmon")) +
  labs(
    title = "SE4 across IGT",
    x = "Trial",
    y = "Frequency",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 0.45), xlim = c(0, 100))

# ---- Combine plots ----

p_combined <- grid.arrange(p_SE3, p_SE4, ncol = 1)

# Save
ggsave("VSE_Sequential_Exploration.png", p_combined, width = 10, height = 12, dpi = 300)

cat("Saved: VSE_Sequential_Exploration.png\n")


# ---- Print summary statistics ----
cat("\n=== SEQUENTIAL EXPLORATION SUMMARY ===\n")
cat("\nSE3 (overall mean across trials):\n")
cat("  Controls:    ", round(mean(SE3_ctr_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Heroin:      ", round(mean(SE3_opi_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Amphetamine: ", round(mean(SE3_amp_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Chance:       0.333\n")

cat("\nSE4 (overall mean across trials):\n")
cat("  Controls:    ", round(mean(SE4_ctr_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Heroin:      ", round(mean(SE4_opi_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Amphetamine: ", round(mean(SE4_amp_sum$mean, na.rm=TRUE), 3), "\n")
cat("  Chance:       0.090\n")

# Create comparison plots
png("VSE_group_comparison_parameters.png", width=1200, height=800)
par(mfrow=c(2,3))

compare_plot("theta (value sensitivity)", 
             Y$mu_theta_ctr, Y$mu_theta_opi, Y$mu_theta_amp, c(0, 1))
compare_plot("delta (decay)", 
             Y$mu_delta_ctr, Y$mu_delta_opi, Y$mu_delta_amp, c(0, 1))
compare_plot("alpha (exploration learning)", 
             Y$mu_alpha_ctr, Y$mu_alpha_opi, Y$mu_alpha_amp, c(0, 1))
compare_plot("phi (exploration bonus)", 
             Y$mu_phi_ctr, Y$mu_phi_opi, Y$mu_phi_amp, c(-2, 2))
compare_plot("c (consistency)", 
             Y$mu_c_ctr, Y$mu_c_opi, Y$mu_c_amp, c(0, 5))

dev.off()
cat("\nSaved: VSE_group_comparison_parameters.png\n")


# ============================================================
# Post-hoc Posterior Predictive Check for VSE
# ============================================================

# ---- Function to simulate VSE model ----
simulate_VSE <- function(theta, delta, alpha, phi, c, x_actual, X_actual, ntrials) {
  # Returns predicted choices for each trial based on actual history
  
  # Initialize
  exploit <- matrix(0, nrow = ntrials, ncol = 4)
  explore <- matrix(0, nrow = ntrials, ncol = 4)
  Ev <- matrix(0, nrow = ntrials, ncol = 4)
  p <- matrix(0.25, nrow = ntrials, ncol = 4)  # t=1 is uniform
  predicted_choices <- rep(NA, ntrials)
  
  # Trial 1 is random (no prediction)
  predicted_choices[1] <- NA
  
  # Simulate trial by trial using ACTUAL choices
  for (t in 2:ntrials) {
    
    # Get previous choice and outcome FROM ACTUAL DATA
    prev_choice <- x_actual[t-1]
    prev_outcome <- X_actual[t-1]
    
    # Skip if previous choice is NA
    if (is.na(prev_choice) || is.na(prev_outcome)) {
      next
    }
    
    # Separate reward and loss
    R <- max(0, prev_outcome)
    L <- abs(min(0, prev_outcome))
    
    # Value function
    v <- R^theta - L^theta
    
    for (d in 1:4) {
      # Exploitation
      if (d == prev_choice) {
        exploit[t, d] <- delta * exploit[t-1, d] + v
      } else {
        exploit[t, d] <- delta * exploit[t-1, d]
      }
      
      # Exploration
      if (d == prev_choice) {
        explore[t, d] <- 0
      } else {
        explore[t, d] <- explore[t-1, d] + alpha * (phi - explore[t-1, d])
      }
      
      # Combined expected value
      Ev[t, d] <- exploit[t, d] + explore[t, d]
    }
    
    # Softmax
    exp_p <- exp(c * Ev[t, ])
    p[t, ] <- exp_p / sum(exp_p)
    
    # Store predicted choice (deterministic: argmax)
    predicted_choices[t] <- which.max(p[t, ])
  }
  
  return(predicted_choices)
}

# ---- Function to compute accuracy for one subject across posterior samples ----
compute_subject_accuracy <- function(theta_samples, delta_samples, alpha_samples, 
                                     phi_samples, c_samples, 
                                     x_true, X_true, ntrials) {
  
  n_samples <- length(theta_samples)
  accuracies <- rep(NA, n_samples)
  
  for (i in 1:n_samples) {
    # Simulate with this posterior sample's parameters
    predicted <- simulate_VSE(
      theta = theta_samples[i],
      delta = delta_samples[i],
      alpha = alpha_samples[i],
      phi = phi_samples[i],
      c = c_samples[i],
      x_actual = x_true,
      X_actual = X_true,
      ntrials = ntrials
    )
    
    # Compare predicted vs actual choices (excluding trial 1)
    correct <- predicted[2:ntrials] == x_true[2:ntrials]
    accuracies[i] <- mean(correct, na.rm = TRUE)
  }
  
  # Return mean accuracy across posterior samples
  return(mean(accuracies))
}

# ---- Apply to all subjects in each group ----

cat("\n=== Computing Post-hoc Posterior Predictive Checks ===\n")

# Controls
cat("Computing for Controls...\n")
accuracy_ctr <- rep(NA, nsubs_ctr)
for (s in 1:nsubs_ctr) {
  accuracy_ctr[s] <- compute_subject_accuracy(
    theta_samples = Y$theta_ctr[, s],
    delta_samples = Y$delta_ctr[, s],
    alpha_samples = Y$alpha_ctr[, s],
    phi_samples = Y$phi_ctr[, s],
    c_samples = Y$c_ctr[, s],
    x_true = x_ctr[s, ],
    X_true = X_ctr[s, ],
    ntrials = ntrials_ctr[s]
  )
  if (s %% 10 == 0) cat("  Processed", s, "/", nsubs_ctr, "subjects\n")
}

# Heroin
cat("Computing for Heroin...\n")
accuracy_opi <- rep(NA, nsubs_opi)
for (s in 1:nsubs_opi) {
  accuracy_opi[s] <- compute_subject_accuracy(
    theta_samples = Y$theta_opi[, s],
    delta_samples = Y$delta_opi[, s],
    alpha_samples = Y$alpha_opi[, s],
    phi_samples = Y$phi_opi[, s],
    c_samples = Y$c_opi[, s],
    x_true = x_opi[s, ],
    X_true = X_opi[s, ],
    ntrials = ntrials_opi[s]
  )
  if (s %% 10 == 0) cat("  Processed", s, "/", nsubs_opi, "subjects\n")
}

# Amphetamine
cat("Computing for Amphetamine...\n")
accuracy_amp <- rep(NA, nsubs_amp)
for (s in 1:nsubs_amp) {
  accuracy_amp[s] <- compute_subject_accuracy(
    theta_samples = Y$theta_amp[, s],
    delta_samples = Y$delta_amp[, s],
    alpha_samples = Y$alpha_amp[, s],
    phi_samples = Y$phi_amp[, s],
    c_samples = Y$c_amp[, s],
    x_true = x_amp[s, ],
    X_true = X_amp[s, ],
    ntrials = ntrials_amp[s]
  )
  if (s %% 10 == 0) cat("  Processed", s, "/", nsubs_amp, "subjects\n")
}

# ---- Print summary ----
cat("\n=== POSTERIOR PREDICTIVE CHECK RESULTS ===\n")
cat("Mean Accuracy (excluding trial 1):\n")
cat("  Controls:    ", round(mean(accuracy_ctr), 3), "\n")
cat("  Heroin:      ", round(mean(accuracy_opi), 3), "\n")
cat("  Amphetamine: ", round(mean(accuracy_amp), 3), "\n")
cat("\nChance level: 0.25\n")

# ---- Plot results ----
png("VSE_posterior_predictive_check.png", width = 1400, height = 500)
par(mfrow = c(1, 3))

barplot(accuracy_ctr * 100, ylim = c(0, 80), 
        main = "Controls (HC)", xlab = "Participants", 
        ylab = "Percent correct", col = "lightgreen", border = NA)
abline(h = 25, lty = 2, lwd = 2)
abline(h = mean(accuracy_ctr) * 100, lwd = 2)
legend("topright", legend = c("Mean accuracy", "Chance level"), 
       lty = c(1, 2), lwd = 2, bty = "n")

barplot(accuracy_opi * 100, ylim = c(0, 80), 
        main = "Heroin (OPI)", xlab = "Participants", 
        ylab = "Percent correct", col = "lightblue", border = NA)
abline(h = 25, lty = 2, lwd = 2)
abline(h = mean(accuracy_opi) * 100, lwd = 2)

barplot(accuracy_amp * 100, ylim = c(0, 80), 
        main = "Amphetamine (AMP)", xlab = "Participants", 
        ylab = "Percent correct", col = "salmon", border = NA)
abline(h = 25, lty = 2, lwd = 2)
abline(h = mean(accuracy_amp) * 100, lwd = 2)

dev.off()
cat("\nSaved: VSE_posterior_predictive_check.png\n")

# ---- Additional diagnostic: Distribution of accuracies ----
png("VSE_PPC_accuracy_distributions.png", width = 1200, height = 400)
par(mfrow = c(1, 3))

hist(accuracy_ctr * 100, breaks = 20, xlim = c(0, 100),
     main = "Controls", xlab = "Accuracy (%)", 
     col = "lightgreen", border = "white")
abline(v = mean(accuracy_ctr) * 100, lwd = 2, col = "darkgreen")
abline(v = 25, lwd = 2, lty = 2, col = "black")

hist(accuracy_opi * 100, breaks = 20, xlim = c(0, 100),
     main = "Heroin", xlab = "Accuracy (%)", 
     col = "lightblue", border = "white")
abline(v = mean(accuracy_opi) * 100, lwd = 2, col = "blue")
abline(v = 25, lwd = 2, lty = 2, col = "black")

hist(accuracy_amp * 100, breaks = 20, xlim = c(0, 100),
     main = "Amphetamine", xlab = "Accuracy (%)", 
     col = "salmon", border = "white")
abline(v = mean(accuracy_amp) * 100, lwd = 2, col = "red")
abline(v = 25, lwd = 2, lty = 2, col = "black")

dev.off()
cat("Saved: VSE_PPC_accuracy_distributions.png\n")

cat("\n=== PPC COMPLETE ===\n")
# ============================================================
# Improved Statistical Testing for 3 Groups
# ============================================================


# ---- Prepare data (assuming you already have the variables) ----
calculate_net_outcome <- function(x_mat, ntrials_vec, nsubs) {
  # rewards_mat is the actual outcome per trial (can be negative for losses)
  net_outcomes <- rep(NA, nsubs)
  
  for (s in 1:nsubs) {
    # take only the valid trials
    outcomes <- x_mat[s, 1:ntrials_vec[s]]
    outcomes <- na.omit(outcomes)
    
    # net outcome = sum of all trial outcomes
    net_outcomes[s] <- sum(outcomes)
  }
  
  return(net_outcomes)
}


net_score_ctr <- calculate_net_outcome(X_ctr, ntrials_ctr, nsubs_ctr)
net_score_opi <- calculate_net_outcome(X_opi, ntrials_opi, nsubs_opi)
net_score_amp <- calculate_net_outcome(X_amp, ntrials_amp, nsubs_amp)



# ---- Calculate overall SE3/SE4 per participant ----

calculate_overall_SE <- function(SE_mat, ntrials_vec, nsubs) {
  overall_SE <- rep(NA, nsubs)
  
  for (s in 1:nsubs) {
    se_values <- SE_mat[s, 1:ntrials_vec[s]]
    overall_SE[s] <- mean(se_values, na.rm = TRUE)
  }
  
  return(overall_SE)
}


# Using the SE matrices we already computed
SE3_overall_ctr <- calculate_overall_SE(SE_ctr$SE3, ntrials_ctr, nsubs_ctr)
SE3_overall_opi <- calculate_overall_SE(SE_opi$SE3, ntrials_opi, nsubs_opi)
SE3_overall_amp <- calculate_overall_SE(SE_amp$SE3, ntrials_amp, nsubs_amp)

SE4_overall_ctr <- calculate_overall_SE(SE_ctr$SE4, ntrials_ctr, nsubs_ctr)
SE4_overall_opi <- calculate_overall_SE(SE_opi$SE4, ntrials_opi, nsubs_opi)
SE4_overall_amp <- calculate_overall_SE(SE_amp$SE4, ntrials_amp, nsubs_amp)


# Create long-format dataframes
df_net_score <- data.frame(
  value = c(net_score_ctr, net_score_opi, net_score_amp),
  group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                     c(nsubs_ctr, nsubs_opi, nsubs_amp)),
                 levels = c("Heroin", "Controls", "Amphetamine"))
)

df_SE3 <- data.frame(
  value = c(SE3_overall_ctr, SE3_overall_opi, SE3_overall_amp),
  group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                     c(nsubs_ctr, nsubs_opi, nsubs_amp)),
                 levels = c("Heroin", "Controls", "Amphetamine"))
)

df_SE4 <- data.frame(
  value = c(SE4_overall_ctr, SE4_overall_opi, SE4_overall_amp),
  group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                     c(nsubs_ctr, nsubs_opi, nsubs_amp)),
                 levels = c("Heroin", "Controls", "Amphetamine"))
)

# ---- Run ANOVA and Post-hoc Tests ----

# Function to run complete statistical analysis
run_stats <- function(df, measure_name) {
  
  cat("\n===============================================\n")
  cat(measure_name, "\n")
  cat("===============================================\n")
  
  # 1. Descriptive statistics
  desc_stats <- df %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = 'drop'
    )
  
  print(desc_stats)
  
  # 2. One-way ANOVA
  anova_result <- aov(value ~ group, data = df)
  anova_summary <- summary(anova_result)
  
  cat("\nOne-way ANOVA:\n")
  print(anova_summary)
  
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]
  
  cat("\nF(", anova_summary[[1]]$Df[1], ", ", anova_summary[[1]]$Df[2], 
      ") = ", round(f_stat, 3), ", p = ", round(p_value, 4), "\n", sep = "")
  
  # 3. Post-hoc tests (only if ANOVA is significant)
  if (p_value < 0.05) {
    cat("\n*** ANOVA significant - Running post-hoc tests ***\n")
    
    # Tukey's HSD (honestly significant difference)
    posthoc <- TukeyHSD(anova_result)
    cat("\nTukey's HSD Post-hoc Comparisons:\n")
    print(posthoc)
    
    # Extract significant comparisons
    posthoc_df <- as.data.frame(posthoc$group)
    posthoc_df$comparison <- rownames(posthoc_df)
    posthoc_df$significant <- posthoc_df$`p adj` < 0.05
    
    return(list(
      anova_p = p_value,
      posthoc = posthoc_df,
      significant = TRUE
    ))
    
  } else {
    cat("\n*** ANOVA not significant - No post-hoc tests needed ***\n")
    
    return(list(
      anova_p = p_value,
      posthoc = NULL,
      significant = FALSE
    ))
  }
}

# Run statistics for each measure
stats_net <- run_stats(df_net_score, "NET SCORE")
stats_SE3 <- run_stats(df_SE3, "SEQUENTIAL EXPLORATION 3")
stats_SE4 <- run_stats(df_SE4, "SEQUENTIAL EXPLORATION 4")

# ---- Enhanced plotting function with proper significance ----

plot_with_significance <- function(df, stats_result, y_label, title, y_limits = NULL) {
  
  # Calculate summary statistics
  summary_stats <- df %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n()),
      .groups = 'drop'
    )
  
  # Define colors
  colors <- c("Heroin" = "lightblue", "Controls" = "lightgreen", "Amphetamine" = "salmon")
  
  p <- ggplot(df, aes(x = group, y = value, color = group, fill = group)) +
    # Bar for mean
    geom_bar(data = summary_stats, aes(x = group, y = mean), 
             stat = "identity", alpha = 0.6, width = 0.7) +
    # Error bars (SD)
    geom_errorbar(data = summary_stats, 
                  aes(x = group, y = mean, ymin = mean - sd, ymax = mean + sd),
                  width = 0.2, color = "black", size = 0.8) +
    # Individual points (jittered)
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.6) +
    # Styling
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      title = title,
      x = "",
      y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 11)
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # Add significance brackets if ANOVA was significant
  if (stats_result$significant && !is.null(stats_result$posthoc)) {
    
    y_max <- max(df$value, na.rm = TRUE)
    y_min <- min(df$value, na.rm = TRUE)
    y_range <- y_max - y_min
    
    # Filter for significant comparisons only
    sig_comparisons <- stats_result$posthoc[stats_result$posthoc$significant, ]
    
    if (nrow(sig_comparisons) > 0) {
      
      y_start <- y_max + y_range * 0.05
      
      for (i in 1:nrow(sig_comparisons)) {
        # Parse comparison string (e.g., "Controls-Heroin")
        comp_name <- sig_comparisons$comparison[i]
        groups <- strsplit(comp_name, "-")[[1]]
        
        x1 <- match(groups[1], levels(df$group))
        x2 <- match(groups[2], levels(df$group))
        
        y_pos <- y_start + (i - 1) * y_range * 0.08
        
        p <- p + 
          annotate("segment", x = x1, xend = x2, 
                   y = y_pos, yend = y_pos, size = 0.5) +
          annotate("segment", x = x1, xend = x1, 
                   y = y_pos - y_range * 0.02, yend = y_pos, size = 0.5) +
          annotate("segment", x = x2, xend = x2, 
                   y = y_pos - y_range * 0.02, yend = y_pos, size = 0.5) +
          annotate("text", x = (x1 + x2) / 2, y = y_pos + y_range * 0.02,
                   label = "*", size = 6)
      }
    }
  }
  
  return(p)
}

# ---- Create plots ----

p1 <- plot_with_significance(df_net_score, stats_net, "Net Score", "Net Score")
p2 <- plot_with_significance(df_SE3, stats_SE3, "Frequency", 
                             "Sequential Exploration 3", y_limits = c(0, 0.6))
p3 <- plot_with_significance(df_SE4, stats_SE4, "Frequency", 
                             "Sequential Exploration 4", y_limits = c(0, 0.4))

# ---- Combine and save ----

p_combined <- grid.arrange(p1, p2, p3, ncol = 3)

ggsave("VSE_Group_Comparisons.png", p_combined, width = 14, height = 5, dpi = 300)

cat("\n\n=== PLOTS SAVED ===\n")
cat("File: VSE_Group_Comparisons.png\n")


#----------Save results----
save(samples, file="VSE_group_comparison_results.RData")
cat("\nSaved: VSE_group_comparison_results.RData\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

