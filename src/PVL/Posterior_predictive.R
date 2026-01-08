# PVL Group Comparison - Posterior Predictive Check
pacman::p_load(ggplot2)

set.seed(1983)

setwd('/work/JoMat/DecisionMaking/')

# --- directories ---
output_dir <- "src/PVL/outputs/posterior_predictive"
fig_dir <- file.path(output_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --- function to calculate maximum of posterior density ---
MPD <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2) return(NA)
  dens <- density(x)
  dens$x[which.max(dens$y)]
}

# --- Load group comparison results ---
cat("\n========================================\n")
cat("Loading PVL Group Comparison Results\n")
cat("========================================\n")

load("src/PVL/outputs/group_comparison/PVL_group_comparison_results.RData")
Y <- samples$BUGSoutput$sims.list

cat("Loaded posterior samples\n")

# --- Define groups and their data files ---
groups <- list(
  control = list(
    name = "Healthy Control",
    data_file = "data/IGTdata_healthy_control_clean.txt",
    p_array = Y$p_ctr,
    color = "lightgreen"
  ),
  heroin = list(
    name = "Heroin",
    data_file = "data/IGTdata_heroin_clean.txt",
    p_array = Y$p_her,
    color = "lightblue"
  ),
  amphetamine = list(
    name = "Amphetamine",
    data_file = "data/IGTdata_amphetamine_clean.txt",
    p_array = Y$p_amp,
    color = "salmon"
  )
)

# --- Loop over groups ---
for (grp in names(groups)) {
  
  cat("\n========================================\n")
  cat("Processing group:", groups[[grp]]$name, "\n")
  cat("========================================\n")
  
  # --- Extract group-specific data ---
  p_full <- groups[[grp]]$p_array
  
  # --- DIAGNOSTIC: Check structure of p array ---
  cat("\n--- Checking posterior structure ---\n")
  cat("Dimensions of p array:", dim(p_full), "\n")
  cat("  [1] = iterations:", dim(p_full)[1], "\n")
  cat("  [2] = subjects:", dim(p_full)[2], "\n")
  cat("  [3] = trials (2:100):", dim(p_full)[3], "\n")
  cat("  [4] = decks:", dim(p_full)[4], "\n")
  cat("Note: Trial 1 cannot be predicted (no prior info), so posterior contains trials 2-100\n")
  
  # --- Load the original data ---
  group_data <- read.table(groups[[grp]]$data_file, header = TRUE)
  
  subIDs <- unique(group_data$subjID)
  nsubs <- length(subIDs)
  ntrials_max <- 100
  
  cat("\nNumber of subjects:", nsubs, "\n")
  
  # Prepare arrays
  x_all <- array(0, c(nsubs, ntrials_max))
  X_all <- array(0, c(nsubs, ntrials_max))
  ntrials_all <- numeric(nsubs)
  
  # --- Loop to pad trials to ntrials_max ---
  x_raw <- group_data$deck
  X_raw <- group_data$gain + group_data$loss
  
  for (s in 1:nsubs) {
    ntrials_all[s] <- length(x_raw[group_data$subjID == subIDs[s]])
    
    x_sub <- x_raw[group_data$subjID == subIDs[s]]
    length(x_sub) <- ntrials_max
    
    X_sub <- X_raw[group_data$subjID == subIDs[s]]
    length(X_sub) <- ntrials_max
    
    x_all[s, ] <- x_sub
    X_all[s, ] <- X_sub
  }
  
  # --- Posterior predictive per subject ---
  pred_success <- numeric(nsubs)
  start_time <- Sys.time()
  
  cat("\n--- Starting posterior predictive checks ---\n")
  
  for (s in 1:nsubs) {
    
    ntrials <- ntrials_all[s]
    if (ntrials == 0) next
    
    x <- x_all[s, 1:ntrials]
    
    # Diagnostic for first subject
    if (s == 1) {
      cat("\n--- Subject 1 diagnostics ---\n")
      cat("Number of trials in data:", ntrials, "\n")
      cat("Number of trials in posterior:", dim(p_full)[3], "\n")
      cat("Note: Posterior has 99 trials (trials 2-100) since trial 1 cannot be predicted\n")
    }
    
    # Check if subject index is valid
    if (s > dim(p_full)[2]) {
      cat("ERROR: Subject", s, "exceeds number of subjects in posterior (", dim(p_full)[2], ")\n")
      next
    }
    
    # The posterior contains trials 2:100 (99 trials total)
    # We'll use all available trials from the posterior
    n_posterior_trials <- dim(p_full)[3]
    
    # Extract subject-level posterior for all available trials
    # p_post dimensions: [iterations, trials (2:100), decks]
    p_post <- p_full[, s, , ]
    
    if (s == 1) {
      cat("Extracted p_post dimensions:", dim(p_post), "\n")
      cat("  [1] = iterations:", dim(p_post)[1], "\n")
      cat("  [2] = trials (2:100):", dim(p_post)[2], "\n")
      cat("  [3] = decks:", dim(p_post)[3], "\n")
    }
    
    # Posterior predictive choices
    # p_post contains predictions for trials 2:100 (indexed as 1:99 in the array)
    x_predict <- numeric(ntrials)
    x_predict[1] <- NA  # No prediction for trial 1
    
    # Loop through trials 2 to 100
    # For trial t, use posterior index t-1 (since posterior starts at trial 2)
    for (t in 2:min(ntrials, n_posterior_trials + 1)) {
      posterior_idx <- t - 1  # Trial 2 -> index 1, Trial 3 -> index 2, etc.
      
      p_predict <- c(
        MPD(p_post[, posterior_idx, 1]),
        MPD(p_post[, posterior_idx, 2]),
        MPD(p_post[, posterior_idx, 3]),
        MPD(p_post[, posterior_idx, 4])
      )
      
      # Normalize to ensure probabilities sum to 1
      p_predict <- p_predict / sum(p_predict)
      
      # Debug output for first subject, first few trials
      if (s == 1 && t <= 5) {
        cat(paste0("Trial ", t, " (posterior index ", posterior_idx, 
                   "): p_predict = ", paste(round(p_predict, 3), collapse = ", "),
                   " | sum = ", round(sum(p_predict), 3), "\n"))
      }
      
      x_predict[t] <- which.max(p_predict)
    }
    
    # Debug output for first subject
    if (s == 1) {
      cat("\nSubject 1 final trial probabilities (trial 100, posterior index 99):\n")
      cat("Mean across decks:\n")
      print(colMeans(p_post[, 99, ]))
      cat("MPD across decks:\n")
      print(apply(p_post[, 99, ], 2, MPD))
      cat("\nFirst 10 actual choices:", x[1:10], "\n")
      cat("First 10 predicted choices:", x_predict[1:10], "\n")
    }
    
    # Count matches (excluding trial 1, comparing trials 2-100)
    # Only compare up to the number of trials we have predictions for
    n_compare_trials <- min(ntrials, n_posterior_trials + 1)
    pred_success[s] <- sum(x_predict[2:n_compare_trials] == x[2:n_compare_trials], na.rm = TRUE)
    
    if (s == 1 || s %% 10 == 0) {
      cat(paste0("Subject ", s, "/", nsubs, " completed (accuracy: ", 
                 round(pred_success[s]/(n_compare_trials - 1), 3), ")\n"))
    }
  }
  
  end_time <- Sys.time()
  cat("\nProcessing time:", difftime(end_time, start_time, units = "secs"), "seconds\n")
  
  # Posterior predictive accuracy
  # We predict trials 2-100 (99 trials), so denominator is 99 (or n_posterior_trials)
  n_posterior_trials <- dim(p_full)[3]
  pred_success_adjust <- pred_success / n_posterior_trials
  
  cat("\n--- Prediction accuracy summary ---\n")
  cat("Mean:", round(mean(pred_success_adjust, na.rm = TRUE), 3), "\n")
  cat("SD:  ", round(sd(pred_success_adjust, na.rm = TRUE), 3), "\n")
  cat("Min: ", round(min(pred_success_adjust, na.rm = TRUE), 3), "\n")
  cat("Max: ", round(max(pred_success_adjust, na.rm = TRUE), 3), "\n")
  
  # Plotting
  pred_df <- data.frame(pred_success_adjust)
  pred_df$sub <- 1:length(pred_success_adjust)
  pred_df$avg <- mean(pred_df$pred_success_adjust, na.rm = TRUE)
  pred_df$std <- sd(pred_df$pred_success_adjust, na.rm = TRUE)
  pred_df$chance <- 0.25
  
  # Calculate accuracy statistics for subtitle
  mean_acc <- mean(pred_df$pred_success_adjust, na.rm = TRUE)
  sd_acc <- sd(pred_df$pred_success_adjust, na.rm = TRUE)
  subtitle_text <- sprintf("Mean Accuracy: %.3f (SD: %.3f)", mean_acc, sd_acc)
  
  p <- ggplot(pred_df, aes(sub, pred_success_adjust)) +
    geom_point(color = groups[[grp]]$color, size = 2.5, alpha = 0.6) +
    geom_line(aes(y = chance), linetype = "dashed", color = "black", size = 0.8) +
    geom_ribbon(aes(xmin = -Inf, xmax = Inf, ymin = avg - std, ymax = avg + std),
                fill = groups[[grp]]$color, alpha = 0.3) +
    geom_line(aes(y = avg), size = 1) +
    ylim(0, 0.8) +
    labs(
      title = paste("Posterior Predictive Accuracy:", groups[[grp]]$name),
      subtitle = subtitle_text,
      x = "Subject",
      y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  plot_path <- file.path(fig_dir, paste0("PPC_", grp, ".png"))
  ggsave(plot_path, plot = p, width = 12, height = 6, dpi = 300)
  cat("\nPlot saved to:", plot_path, "\n")
  
  # Save results
  results_path <- file.path(output_dir, paste0("PPC_results_", grp, ".RData"))
  save(pred_success_adjust, pred_df, file = results_path)
  cat("Results saved to:", results_path, "\n")
  
  cat(paste0("\nFinished group: ", grp, "\n"))
}

# --- Create combined comparison plot ---
cat("\n========================================\n")
cat("Creating Combined Comparison Plot\n")
cat("========================================\n")

# Load all results
results_list <- list()
for (grp in names(groups)) {
  results_path <- file.path(output_dir, paste0("PPC_results_", grp, ".RData"))
  load(results_path)
  
  results_list[[grp]] <- data.frame(
    accuracy = pred_success_adjust,
    group = groups[[grp]]$name,
    color = groups[[grp]]$color
  )
}

# Combine into single dataframe
combined_df <- do.call(rbind, results_list)
combined_df$group <- factor(combined_df$group, 
                            levels = c("Heroin", "Healthy Control", "Amphetamine"))

# Create combined plot
p_combined <- ggplot(combined_df, aes(x = group, y = accuracy, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Heroin" = "lightblue", 
                               "Healthy Control" = "lightgreen", 
                               "Amphetamine" = "salmon")) +
  labs(
    title = "Posterior Predictive Accuracy Comparison",
    subtitle = "Dashed line = chance level (0.25)",
    x = "",
    y = "Proportion Correctly Predicted"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  ylim(0, 0.8)

combined_path <- file.path(fig_dir, "PPC_combined_comparison.png")
ggsave(combined_path, plot = p_combined, width = 10, height = 6, dpi = 300)
cat("Combined plot saved to:", combined_path, "\n")

# Print summary statistics
cat("\n--- Summary Statistics ---\n")
for (grp in names(groups)) {
  grp_data <- combined_df[combined_df$group == groups[[grp]]$name, ]
  cat("\n", groups[[grp]]$name, ":\n", sep = "")
  cat("  Mean:   ", round(mean(grp_data$accuracy, na.rm = TRUE), 3), "\n")
  cat("  SD:     ", round(sd(grp_data$accuracy, na.rm = TRUE), 3), "\n")
  cat("  Median: ", round(median(grp_data$accuracy, na.rm = TRUE), 3), "\n")
  cat("  Range:  ", round(min(grp_data$accuracy, na.rm = TRUE), 3), 
      "-", round(max(grp_data$accuracy, na.rm = TRUE), 3), "\n")
}

cat("\n========================================\n")
cat("✅ ALL GROUPS PROCESSED\n")
cat("========================================\n")