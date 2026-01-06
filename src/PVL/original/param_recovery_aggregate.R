#!/usr/bin/env Rscript

# Aggregates results from all parallel parameter recovery batches

pacman::p_load(ggpubr)

setwd('/work/JoMat/DecisionMaking/src/PVL/original')
source('recov_plot.R')

cat("Aggregating parameter recovery results...\n")

# Initialize aggregation arrays
niterations_total <- 0
all_true_mu_w <- c()
all_true_mu_A <- c()
all_true_mu_theta <- c()
all_true_mu_a <- c()

all_infer_mu_w <- c()
all_infer_mu_A <- c()
all_infer_mu_theta <- c()
all_infer_mu_a <- c()

all_true_lambda_w <- c()
all_true_lambda_A <- c()
all_true_lambda_theta <- c()
all_true_lambda_a <- c()

all_infer_lambda_w <- c()
all_infer_lambda_A <- c()
all_infer_lambda_theta <- c()
all_infer_lambda_a <- c()

# Load and combine results from all batch files
result_files <- list.files('results/', pattern = 'param_recovery_batch_.*\\.RData$', full.names = TRUE)
result_files <- sort(result_files)

if (length(result_files) == 0) {
  stop("No result files found in results/ directory!")
}

cat("Found", length(result_files), "result files\n")

for (file in result_files) {
  cat("Loading:", file, "\n")
  load(file)
  
  # Combine results
  all_true_mu_w <- c(all_true_mu_w, true_mu_w)
  all_true_mu_A <- c(all_true_mu_A, true_mu_A)
  all_true_mu_theta <- c(all_true_mu_theta, true_mu_theta)
  all_true_mu_a <- c(all_true_mu_a, true_mu_a)
  
  all_infer_mu_w <- c(all_infer_mu_w, infer_mu_w)
  all_infer_mu_A <- c(all_infer_mu_A, infer_mu_A)
  all_infer_mu_theta <- c(all_infer_mu_theta, infer_mu_theta)
  all_infer_mu_a <- c(all_infer_mu_a, infer_mu_a)
  
  all_true_lambda_w <- c(all_true_lambda_w, true_lambda_w)
  all_true_lambda_A <- c(all_true_lambda_A, true_lambda_A)
  all_true_lambda_theta <- c(all_true_lambda_theta, true_lambda_theta)
  all_true_lambda_a <- c(all_true_lambda_a, true_lambda_a)
  
  all_infer_lambda_w <- c(all_infer_lambda_w, infer_lambda_w)
  all_infer_lambda_A <- c(all_infer_lambda_A, infer_lambda_A)
  all_infer_lambda_theta <- c(all_infer_lambda_theta, infer_lambda_theta)
  all_infer_lambda_a <- c(all_infer_lambda_a, infer_lambda_a)
}

cat("Total iterations:", length(all_true_mu_w), "\n")

# Create recovery plots for means
cat("Creating plots for means...\n")
pl1 <- recov_plot(all_true_mu_w, all_infer_mu_w, c("true mu_w", "infer mu_w"), 'smoothed linear fit')
pl2 <- recov_plot(all_true_mu_A, all_infer_mu_A, c("true mu_A", "infer mu_A"), 'smoothed linear fit')
pl3 <- recov_plot(all_true_mu_theta, all_infer_mu_theta, c("true mu_theta", "infer mu_theta"), 'smoothed linear fit')
pl4 <- recov_plot(all_true_mu_a, all_infer_mu_a, c("true mu_a", "infer mu_a"), 'smoothed linear fit')

png('results/recovery_plots_means.png', width = 1200, height = 1000)
ggarrange(pl1, pl2, pl3, pl4)
dev.off()
cat("Saved: results/recovery_plots_means.png\n")

# Create recovery plots for precision
cat("Creating plots for precision...\n")
pl1 <- recov_plot(all_infer_lambda_w, 1/(all_true_lambda_w^2), c("infer lambda_w", "true lambda_w"), 'smoothed linear fit')
pl2 <- recov_plot(all_infer_lambda_A, 1/(all_true_lambda_A^2), c("infer lambda_A", "true lambda_A"), 'smoothed linear fit')
pl3 <- recov_plot(all_infer_lambda_theta, 1/(all_true_lambda_theta^2), c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl4 <- recov_plot(all_infer_lambda_a, 1/(all_true_lambda_a^2), c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')

png('results/recovery_plots_precision.png', width = 1200, height = 1000)
ggarrange(pl1, pl2, pl3, pl4)
dev.off()
cat("Saved: results/recovery_plots_precision.png\n")

# Create recovery plots for SD
cat("Creating plots for SD...\n")
pl1 <- recov_plot(1/sqrt(all_infer_lambda_w), all_true_lambda_w, c("infer lambda_w", "true lambda_w"), 'smoothed linear fit')
pl2 <- recov_plot(1/sqrt(all_infer_lambda_A), all_true_lambda_A, c("infer lambda_A", "true lambda_A"), 'smoothed linear fit')
pl3 <- recov_plot(1/sqrt(all_infer_lambda_theta), all_true_lambda_theta, c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl4 <- recov_plot(1/sqrt(all_infer_lambda_a), all_true_lambda_a, c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')

png('results/recovery_plots_sd.png', width = 1200, height = 1000)
ggarrange(pl1, pl2, pl3, pl4)
dev.off()
cat("Saved: results/recovery_plots_sd.png\n")

# Save aggregated results
save(all_true_mu_w, all_true_mu_A, all_true_mu_theta, all_true_mu_a,
     all_infer_mu_w, all_infer_mu_A, all_infer_mu_theta, all_infer_mu_a,
     all_true_lambda_w, all_true_lambda_A, all_true_lambda_theta, all_true_lambda_a,
     all_infer_lambda_w, all_infer_lambda_A, all_infer_lambda_theta, all_infer_lambda_a,
     file = 'results/aggregated_recovery_results.RData')
cat("Saved: results/aggregated_recovery_results.RData\n")

cat("\nAggregation complete!\n")