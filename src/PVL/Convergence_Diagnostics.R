library(R2jags)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

setwd('/work/JoMat/DecisionMaking/')

# Create output folder
output_dir <- "src/PVL/outputs/group_comparison/convergence_diagnostics"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to create trace plot for a parameter
create_trace_plot <- function(param_array, param_name, rhat_val, chains = 3) {
  n_iter <- nrow(param_array)
  
  trace_df <- data.frame(
    iteration = rep(1:n_iter, chains),
    value = as.vector(param_array),
    chain = factor(rep(1:chains, each = n_iter))
  )
  
  p <- ggplot(trace_df, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    labs(title = paste("Trace Plot:", param_name),
         subtitle = paste("R-hat =", round(rhat_val, 4)),
         x = "Iteration", y = "Parameter Value") +
    scale_color_manual(values = c("1" = "blue", "2" = "red", "3" = "green")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Load group comparison results
cat("\n========================================\n")
cat("Loading PVL Group Comparison Results\n")
cat("========================================\n\n")

load("src/PVL/outputs/group_comparison/PVL_group_comparison_results.RData")

# Extract MCMC samples and R-hat values
Y <- samples$BUGSoutput$sims.list
mcmc_array <- samples$BUGSoutput$sims.array
rhat <- samples$BUGSoutput$summary[, "Rhat"]

cat("Loaded posterior samples from group comparison model\n")
cat("Total parameters:", length(rhat), "\n\n")

# =============================
# MODEL SUMMARY & DIC
# =============================

cat("========================================\n")
cat("MODEL SUMMARY\n")
cat("========================================\n\n")

# Extract DIC
dic <- samples$BUGSoutput$DIC
pD <- samples$BUGSoutput$pD
deviance <- samples$BUGSoutput$mean$deviance

cat("Deviance Information Criterion (DIC):", round(dic, 2), "\n")
cat("Effective number of parameters (pD):", round(pD, 2), "\n")
cat("Mean deviance:", round(deviance, 2), "\n\n")

cat("MCMC Settings:\n")
cat("  Chains:", samples$BUGSoutput$n.chains, "\n")
cat("  Iterations per chain:", samples$BUGSoutput$n.iter, "\n")
cat("  Burn-in:", samples$BUGSoutput$n.burnin, "\n")
cat("  Thin:", samples$BUGSoutput$n.thin, "\n")
cat("  Total samples:", samples$BUGSoutput$n.sims, "\n\n")

# Save model summary
model_summary_file <- file.path(output_dir, "model_summary.txt")
sink(model_summary_file)

cat("=======================================================\n")
cat("         PVL GROUP COMPARISON MODEL SUMMARY          \n")
cat("=======================================================\n")
cat("\nDate:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== MODEL FIT ===\n")
cat("\nDeviance Information Criterion (DIC):", round(dic, 2), "\n")
cat("Effective number of parameters (pD):", round(pD, 2), "\n")
cat("Mean deviance:", round(deviance, 2), "\n")

cat("\n=== MCMC SETTINGS ===\n")
cat("Chains:", samples$BUGSoutput$n.chains, "\n")
cat("Iterations per chain:", samples$BUGSoutput$n.iter, "\n")
cat("Burn-in:", samples$BUGSoutput$n.burnin, "\n")
cat("Thin:", samples$BUGSoutput$n.thin, "\n")
cat("Total samples:", samples$BUGSoutput$n.sims, "\n")

cat("\n=== DATA SUMMARY ===\n")
# Extract from parameter dimensions
n_ctr <- dim(Y$p_ctr)[2]
n_her <- dim(Y$p_her)[2]
n_amp <- dim(Y$p_amp)[2]

cat("Controls:", n_ctr, "subjects\n")
cat("Heroin:", n_her, "subjects\n")
cat("Amphetamine:", n_amp, "subjects\n")
cat("Total subjects:", n_ctr + n_her + n_amp, "\n")

sink()
cat("Model summary saved to:", model_summary_file, "\n\n")

# =============================
# CONVERGENCE DIAGNOSTICS
# =============================

# Define groups and their parameters
groups <- list(
  control = list(
    name = "Healthy Control",
    params = c("mu_w_ctr", "mu_A_ctr", "mu_theta_ctr", "mu_a_ctr",
               "lambda_w_ctr", "lambda_A_ctr", "lambda_theta_ctr", "lambda_a_ctr"),
    color = "lightgreen"
  ),
  heroin = list(
    name = "Heroin",
    params = c("mu_w_her", "mu_A_her", "mu_theta_her", "mu_a_her",
               "lambda_w_her", "lambda_A_her", "lambda_theta_her", "lambda_a_her"),
    color = "lightblue"
  ),
  amphetamine = list(
    name = "Amphetamine",
    params = c("mu_w_amp", "mu_A_amp", "mu_theta_amp", "mu_a_amp",
               "lambda_w_amp", "lambda_A_amp", "lambda_theta_amp", "lambda_a_amp"),
    color = "salmon"
  )
)

# Calculate overall convergence stats
max_rhat <- max(rhat, na.rm = TRUE)
n_fail <- sum(rhat > 1.1, na.rm = TRUE)
pct_fail <- round(100 * n_fail / length(rhat), 2)

# Loop through each group
for (grp in names(groups)) {
  
  cat("========================================\n")
  cat("Processing:", groups[[grp]]$name, "\n")
  cat("========================================\n\n")
  
  group_params <- groups[[grp]]$params
  
  # =============================
  # 1. Trace Plots - Group Level
  # =============================
  
  # Create trace plots for group-level parameters
  trace_plots <- list()
  
  for (param in group_params) {
    if (param %in% dimnames(mcmc_array)[[3]]) {
      trace_plots[[param]] <- create_trace_plot(mcmc_array[, , param], param, rhat[param])
    }
  }
  
  # Combine trace plots for means
  means_params <- group_params[1:4]
  means_plots <- trace_plots[means_params]
  p_traces_means <- plot_grid(plotlist = means_plots, ncol = 2)
  ggsave(file.path(output_dir, paste0("traces_group_means_", grp, ".png")),
         plot = p_traces_means, width = 14, height = 10, dpi = 300)
  
  # Combine trace plots for precisions
  precision_params <- group_params[5:8]
  precision_plots <- trace_plots[precision_params]
  p_traces_precision <- plot_grid(plotlist = precision_plots, ncol = 2)
  ggsave(file.path(output_dir, paste0("traces_group_precisions_", grp, ".png")),
         plot = p_traces_precision, width = 14, height = 10, dpi = 300)
  
  cat("Trace plots saved for:", groups[[grp]]$name, "\n")
  
  # =============================
  # 2. Save Summary Report
  # =============================
  
  summary_file <- file.path(output_dir, paste0("convergence_summary_", grp, ".txt"))
  sink(summary_file)
  
  cat("=======================================================\n")
  cat("         MCMC CONVERGENCE DIAGNOSTICS REPORT          \n")
  cat("=======================================================\n")
  cat("\nGroup:", groups[[grp]]$name, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  cat("\n=== GROUP-LEVEL PARAMETERS ===\n")
  cat("\nParameter          R-hat\n")
  cat("--------------------------------\n")
  for (param in group_params) {
    if (param %in% names(rhat)) {
      cat(sprintf("%-17s  %.4f", param, rhat[param]))
      if (rhat[param] > 1.1) cat("  ⚠")
      if (rhat[param] > 1.2) cat("  ✗")
      cat("\n")
    }
  }
  
  # Get subject-level parameters for this group
  if (grp == "control") {
    subject_params <- grep("^w_ctr\\[|^A_ctr\\[|^theta_ctr\\[|^a_ctr\\[", names(rhat), value = TRUE)
  } else if (grp == "heroin") {
    subject_params <- grep("^w_her\\[|^A_her\\[|^theta_her\\[|^a_her\\[", names(rhat), value = TRUE)
  } else if (grp == "amphetamine") {
    subject_params <- grep("^w_amp\\[|^A_amp\\[|^theta_amp\\[|^a_amp\\[", names(rhat), value = TRUE)
  }
  
  if (length(subject_params) > 0) {
    subject_rhat <- rhat[subject_params]
    
    cat("\n=== SUBJECT-LEVEL PARAMETERS ===\n")
    cat("\nSummary statistics for subject-level R-hat:\n")
    cat("  Mean:  ", round(mean(subject_rhat, na.rm = TRUE), 4), "\n")
    cat("  Median:", round(median(subject_rhat, na.rm = TRUE), 4), "\n")
    cat("  Min:   ", round(min(subject_rhat, na.rm = TRUE), 4), "\n")
    cat("  Max:   ", round(max(subject_rhat, na.rm = TRUE), 4), "\n")
    cat("  SD:    ", round(sd(subject_rhat, na.rm = TRUE), 4), "\n")
    
    cat("\nN of parameters with R-hat > 1.1:", sum(subject_rhat > 1.1, na.rm = TRUE), "\n")
    cat("N of parameters with R-hat > 1.2:", sum(subject_rhat > 1.2, na.rm = TRUE), "\n")
  }
  
  sink()
  
  cat("Summary saved for:", groups[[grp]]$name, "\n\n")
}

# =============================
# OVERALL CONVERGENCE REPORT
# =============================

cat("========================================\n")
cat("Creating Overall Convergence Report\n")
cat("========================================\n\n")

overall_file <- file.path(output_dir, "convergence_summary_overall.txt")
sink(overall_file)

cat("=======================================================\n")
cat("    OVERALL CONVERGENCE ASSESSMENT - ALL GROUPS       \n")
cat("=======================================================\n")
cat("\nDate:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("\n=== CONVERGENCE CRITERIA ===\n")
cat("Stringent: R-hat < 1.1\n")
cat("Lenient: R-hat < 1.2 (Brooks & Gelman)\n")

cat("\n=== OVERALL STATISTICS ===\n")
cat("Total parameters:", length(rhat), "\n")
cat("Maximum R-hat across all parameters:", round(max_rhat, 4), "\n")
cat("Parameters failing stringent criterion (>1.1):", n_fail, "/", length(rhat), 
    "(", pct_fail, "%)\n")
cat("Parameters failing lenient criterion (>1.2):", sum(rhat > 1.2, na.rm = TRUE), "/", length(rhat), "\n")

if (max_rhat < 1.1) {
  cat("\n✓ CONVERGENCE ACHIEVED\n")
  cat("All parameters meet the stringent R-hat < 1.1 criterion.\n")
} else if (max_rhat < 1.2) {
  cat("\n⚠ PARTIAL CONVERGENCE\n")
  cat("Some parameters exceed 1.1 but all are below 1.2.\n")
  cat("Consider running longer chains if precision is critical.\n")
} else {
  cat("\n✗ CONVERGENCE ISSUES DETECTED\n")
  cat("Some parameters exceed R-hat = 1.2.\n")
  cat("RECOMMENDATION: Increase iterations and/or burnin period.\n")
}

cat("\n=== R-HAT DISTRIBUTION ===\n")
cat("Mean R-hat:", round(mean(rhat, na.rm = TRUE), 4), "\n")
cat("Median R-hat:", round(median(rhat, na.rm = TRUE), 4), "\n")
cat("SD R-hat:", round(sd(rhat, na.rm = TRUE), 4), "\n")

cat("\n=== GROUP DIFFERENCES PARAMETERS ===\n")
diff_params <- grep("^diff_", names(rhat), value = TRUE)
if (length(diff_params) > 0) {
  cat("\nParameter                    R-hat\n")
  cat("-------------------------------------------\n")
  for (param in diff_params) {
    cat(sprintf("%-27s  %.4f", param, rhat[param]))
    if (rhat[param] > 1.1) cat("  ⚠")
    if (rhat[param] > 1.2) cat("  ✗")
    cat("\n")
  }
}

cat("\n=== MODEL FIT ===\n")
cat("DIC:", round(dic, 2), "\n")
cat("pD:", round(pD, 2), "\n")
cat("Mean deviance:", round(deviance, 2), "\n")

sink()

cat("Overall convergence report saved to:", overall_file, "\n")

# =============================
# R-HAT DISTRIBUTION PLOT
# =============================

cat("\n========================================\n")
cat("Creating R-hat Distribution Plots\n")
cat("========================================\n\n")

# Create histogram of all R-hat values
rhat_df <- data.frame(rhat = rhat[!is.na(rhat)])

p_rhat_hist <- ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed", color = "orange", size = 1) +
  geom_vline(xintercept = 1.2, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.1, y = Inf, label = "R-hat = 1.1", 
           vjust = 2, hjust = -0.1, color = "orange") +
  annotate("text", x = 1.2, y = Inf, label = "R-hat = 1.2", 
           vjust = 4, hjust = -0.1, color = "red") +
  labs(title = "Distribution of R-hat Values Across All Parameters",
       subtitle = paste("Max R-hat =", round(max_rhat, 4)),
       x = "R-hat", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "rhat_distribution.png"),
       plot = p_rhat_hist, width = 10, height = 6, dpi = 300)

cat("R-hat distribution plot saved\n")

# =============================
# COMPLETION MESSAGE
# =============================

cat("\n========================================\n")
cat("✅ ALL DIAGNOSTICS COMPLETED\n")
cat("========================================\n")
cat("\nOutputs saved in:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  - model_summary.txt\n")
cat("  - convergence_summary_overall.txt\n")
cat("  - convergence_summary_[group].txt (x3)\n")
cat("  - traces_group_means_[group].png (x3)\n")
cat("  - traces_group_precisions_[group].png (x3)\n")
cat("  - rhat_distribution.png\n\n")