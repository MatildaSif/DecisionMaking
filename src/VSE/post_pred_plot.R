library(R2jags)
library(ggplot2)
library(gridExtra)
library(grid)

set.seed(1702)

setwd('/work/JoMat/DecisionMaking/src/VSE')

cat("\n=== VSE Posterior Predictive Check Analysis ===\n")

#----------Load saved results----
cat("Loading saved model results...\n")
load("VSE_group_comparison_results.RData")

#----------Load original data----
ctr_data <- read.table("../../data/IGTdata_healthy_control_clean.txt", header=TRUE)
opi_data <- read.table("../../data/IGTdata_heroin_clean.txt", header=TRUE)
amp_data <- read.table("../../data/IGTdata_amphetamine_clean.txt", header=TRUE)

#----------Prepare data matrices (same as original script)----
subIDs_ctr <- unique(ctr_data$subjID)
nsubs_ctr <- length(subIDs_ctr)
subIDs_opi <- unique(opi_data$subjID)
nsubs_opi <- length(subIDs_opi)
subIDs_amp <- unique(amp_data$subjID)
nsubs_amp <- length(subIDs_amp)

ntrials_max <- 100

x_raw_ctr <- ctr_data$deck
X_raw_ctr <- (ctr_data$gain + ctr_data$loss)
x_raw_opi <- opi_data$deck
X_raw_opi <- (opi_data$gain + opi_data$loss)
x_raw_amp <- amp_data$deck
X_raw_amp <- (amp_data$gain + amp_data$loss)

# Create matrices
ntrials_ctr <- array(0, c(nsubs_ctr))
x_ctr <- array(0, c(nsubs_ctr, ntrials_max))
X_ctr <- array(0, c(nsubs_ctr, ntrials_max))

ntrials_opi <- array(0, c(nsubs_opi))
x_opi <- array(0, c(nsubs_opi, ntrials_max))
X_opi <- array(0, c(nsubs_opi, ntrials_max))

ntrials_amp <- array(0, c(nsubs_amp))
x_amp <- array(0, c(nsubs_amp, ntrials_max))
X_amp <- array(0, c(nsubs_amp, ntrials_max))

# Fill matrices - Controls
for (s in 1:nsubs_ctr) {
  ntrials_ctr[s] <- length(x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]])
  x_sub_ctr <- x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(x_sub_ctr) <- ntrials_max
  X_sub_ctr <- X_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(X_sub_ctr) <- ntrials_max
  x_ctr[s,] <- x_sub_ctr
  X_ctr[s,] <- X_sub_ctr
}

# Fill matrices - Heroin
for (s in 1:nsubs_opi) {
  ntrials_opi[s] <- length(x_raw_opi[opi_data$subjID == subIDs_opi[s]])
  x_sub_opi <- x_raw_opi[opi_data$subjID == subIDs_opi[s]] 
  length(x_sub_opi) <- ntrials_max
  X_sub_opi <- X_raw_opi[opi_data$subjID == subIDs_opi[s]] 
  length(X_sub_opi) <- ntrials_max
  x_opi[s,] <- x_sub_opi
  X_opi[s,] <- X_sub_opi
}

# Fill matrices - Amphetamine
for (s in 1:nsubs_amp) {
  ntrials_amp[s] <- length(x_raw_amp[amp_data$subjID == subIDs_amp[s]])
  x_sub_amp <- x_raw_amp[amp_data$subjID == subIDs_amp[s]] 
  length(x_sub_amp) <- ntrials_max
  X_sub_amp <- X_raw_amp[amp_data$subjID == subIDs_amp[s]] 
  length(X_sub_amp) <- ntrials_max
  x_amp[s,] <- x_sub_amp
  X_amp[s,] <- X_sub_amp
}

#----------Extract posterior samples----
Y <- samples$BUGSoutput$sims.list

#----------VSE simulation function----
simulate_VSE <- function(theta, delta, alpha, phi, c, x_actual, X_actual, ntrials) {
  exploit <- matrix(0, nrow = ntrials, ncol = 4)
  explore <- matrix(0, nrow = ntrials, ncol = 4)
  Ev <- matrix(0, nrow = ntrials, ncol = 4)
  p <- matrix(0.25, nrow = ntrials, ncol = 4)
  predicted_choices <- rep(NA, ntrials)
  
  predicted_choices[1] <- NA
  
  for (t in 2:ntrials) {
    prev_choice <- x_actual[t-1]
    prev_outcome <- X_actual[t-1]
    
    if (is.na(prev_choice) || is.na(prev_outcome)) {
      next
    }
    
    R <- max(0, prev_outcome)
    L <- abs(min(0, prev_outcome))
    v <- R^theta - L^theta
    
    for (d in 1:4) {
      if (d == prev_choice) {
        exploit[t, d] <- delta * exploit[t-1, d] + v
      } else {
        exploit[t, d] <- delta * exploit[t-1, d]
      }
      
      if (d == prev_choice) {
        explore[t, d] <- 0
      } else {
        explore[t, d] <- explore[t-1, d] + alpha * (phi - explore[t-1, d])
      }
      
      Ev[t, d] <- exploit[t, d] + explore[t, d]
    }
    
    exp_p <- exp(c * Ev[t, ])
    p[t, ] <- exp_p / sum(exp_p)
    predicted_choices[t] <- which.max(p[t, ])
  }
  
  return(predicted_choices)
}

#----------Compute accuracy for one subject----
compute_subject_accuracy <- function(theta_samples, delta_samples, alpha_samples, 
                                     phi_samples, c_samples, 
                                     x_true, X_true, ntrials) {
  n_samples <- length(theta_samples)
  accuracies <- rep(NA, n_samples)
  
  for (i in 1:n_samples) {
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
    
    correct <- predicted[2:ntrials] == x_true[2:ntrials]
    accuracies[i] <- mean(correct, na.rm = TRUE)
  }
  
  return(mean(accuracies))
}

#----------Compute accuracies for all groups----
cat("\nComputing posterior predictive accuracies...\n")

# Controls
cat("  Processing Controls...")
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
}
cat(" Done.\n")

# Heroin
cat("  Processing Heroin...")
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
}
cat(" Done.\n")

# Amphetamine
cat("  Processing Amphetamine...")
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
}
cat(" Done.\n")

#----------Summary statistics----
cat("\n=== POSTERIOR PREDICTIVE CHECK RESULTS ===\n")
cat("\nControls (n =", nsubs_ctr, "):\n")
cat("  Mean accuracy: ", round(mean(accuracy_ctr) * 100, 2), "%\n", sep="")
cat("  SD:            ", round(sd(accuracy_ctr) * 100, 2), "%\n", sep="")
cat("  Range:         ", round(min(accuracy_ctr) * 100, 2), "% - ", 
    round(max(accuracy_ctr) * 100, 2), "%\n", sep="")

cat("\nHeroin (n =", nsubs_opi, "):\n")
cat("  Mean accuracy: ", round(mean(accuracy_opi) * 100, 2), "%\n", sep="")
cat("  SD:            ", round(sd(accuracy_opi) * 100, 2), "%\n", sep="")
cat("  Range:         ", round(min(accuracy_opi) * 100, 2), "% - ", 
    round(max(accuracy_opi) * 100, 2), "%\n", sep="")

cat("\nAmphetamine (n =", nsubs_amp, "):\n")
cat("  Mean accuracy: ", round(mean(accuracy_amp) * 100, 2), "%\n", sep="")
cat("  SD:            ", round(sd(accuracy_amp) * 100, 2), "%\n", sep="")
cat("  Range:         ", round(min(accuracy_amp) * 100, 2), "% - ", 
    round(max(accuracy_amp) * 100, 2), "%\n", sep="")

cat("\nChance level: 25%\n")

#----------Create enhanced barplot with stats box----

create_ppc_barplot <- function(accuracy_vec, group_name, group_color, nsubs) {
  
  # Calculate statistics (keep as proportions, not percentages)
  mean_acc_pct <- mean(accuracy_vec) * 100
  sd_acc_pct <- sd(accuracy_vec) * 100
  mean_acc <- mean(accuracy_vec)
  sd_acc <- sd(accuracy_vec)
  
  # Create barplot with larger text
  bp <- barplot(accuracy_vec * 100, 
                ylim = c(0, 85), 
                main = paste0(group_name, " (n=", nsubs, ")"),
                xlab = "Participants", 
                ylab = "Percent correct",
                col = group_color, 
                border = NA,
                las = 1,
                cex.main = 1.8,
                cex.lab = 1.5,
                cex.axis = 1.3)
  
  # Add reference lines
  abline(h = 25, lty = 2, lwd = 3, col = "gray30")
  abline(h = mean_acc_pct, lwd = 3.5, col = "black")
  
  # Add text box with white background
  rect(xleft = max(bp) * 0.60, 
       xright = max(bp) * 0.98,
       ybottom = 70, 
       ytop = 82,
       col = "white", 
       border = "gray50",
       lwd = 2)
  
  text(x = max(bp) * 0.79, 
       y = 78,
       labels = sprintf("Mean: %.2f", mean_acc),
       cex = 1.6, 
       font = 2)
  
  text(x = max(bp) * 0.79, 
       y = 73.5,
       labels = sprintf("SD: %.2f", sd_acc),
       cex = 1.6, 
       font = 1)
  
  # Add legend
  legend("topleft", 
         legend = c("Mean accuracy", "Chance level"),
         lty = c(1, 2), 
         lwd = c(3.5, 3),
         col = c("black", "gray30"),
         bty = "n",
         cex = 1.3)
}

#----------Generate plots----

# PLOT: Barplots with shaded SD regions (2x2 grid)
png("VSE_PPC_barplots_with_shaded_SD.png", width = 1200, height = 1200)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Function to create barplot with shaded SD
create_ppc_barplot_shaded <- function(accuracy_vec, group_name, group_color, dark_color, nsubs) {
  
  # Calculate statistics
  mean_acc <- mean(accuracy_vec)
  sd_acc <- sd(accuracy_vec)
  mean_pct <- mean_acc * 100
  sd_lower <- (mean_acc - 2*sd_acc) * 100
  sd_upper <- (mean_acc + 2*sd_acc) * 100
  
  # Create barplot
  bp <- barplot(accuracy_vec * 100, 
                ylim = c(0, 85), 
                main = paste0(group_name, " (n=", nsubs, ")"),
                xlab = "Participants", 
                ylab = "Percent correct",
                col = group_color, 
                border = NA,
                las = 1,
                cex.main = 1.8,
                cex.lab = 1.5,
                cex.axis = 1.3)
  
  # Add shaded SD region
  rect(xleft = par("usr")[1], 
       xright = par("usr")[2],
       ybottom = sd_lower, 
       ytop = sd_upper,
       col = adjustcolor(dark_color, alpha.f = 0.2), 
       border = NA)
  
  # Add mean line
  abline(h = mean_pct, lwd = 3.5, col = dark_color)
  
  # Add chance line
  abline(h = 25, lty = 2, lwd = 3, col = "gray30")
  
  # Add text box with stats
  rect(xleft = max(bp) * 0.60, 
       xright = max(bp) * 0.98,
       ybottom = 70, 
       ytop = 82,
       col = "white", 
       border = "gray50",
       lwd = 2)
  
  text(x = max(bp) * 0.79, 
       y = 78,
       labels = sprintf("Mean: %.2f", mean_acc),
       cex = 1.6, 
       font = 2)
  
  text(x = max(bp) * 0.79, 
       y = 73.5,
       labels = sprintf("SD: %.2f", sd_acc),
       cex = 1.6, 
       font = 1)
  
  # Add legend
  legend("topleft", 
         legend = c("Mean", "±1 SD", "Chance"),
         lty = c(1, NA, 2), 
         lwd = c(3.5, NA, 3),
         pch = c(NA, 15, NA),
         col = c(dark_color, adjustcolor(dark_color, alpha.f = 0.3), "gray30"),
         pt.cex = 2,
         bty = "n",
         cex = 1.3)
}

# Create plots
create_ppc_barplot_shaded(accuracy_ctr, "Controls", "lightgreen", "darkgreen", nsubs_ctr)
create_ppc_barplot_shaded(accuracy_opi, "Heroin", "lightblue", "blue", nsubs_opi)
create_ppc_barplot_shaded(accuracy_amp, "Amphetamine", "salmon", "red", nsubs_amp)

# Leave 4th panel empty
plot.new()

dev.off()
cat("Saved: VSE_PPC_barplots_with_shaded_SD.png\n")

#----------Save results to CSV----
results_df <- data.frame(
  Subject_ID = 1:(nsubs_ctr + nsubs_opi + nsubs_amp),
  Group = c(rep("Controls", nsubs_ctr),
            rep("Heroin", nsubs_opi),
            rep("Amphetamine", nsubs_amp)),
  Accuracy = c(accuracy_ctr, accuracy_opi, accuracy_amp) * 100
)

write.csv(results_df, "VSE_PPC_accuracies.csv", row.names = FALSE)
cat("Saved: VSE_PPC_accuracies.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")