# Load the PVL recovery results
setwd('/work/JoMat/DecisionMaking/')

# Create output folder structure
output_dir <- "src/PVL/original/results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load recovery results from output directory
load(file.path(output_dir, "aggregated_recovery_results.RData"))

# Install/load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(corrplot, ggplot2, reshape2)

# Open text file for detailed output
analysis_file <- file.path(output_dir, "recovery_analysis_detailed.txt")
sink(analysis_file, append = FALSE)

# ===== PARAMETER RECOVERY QUALITY (PRIMARY ANALYSIS) =====
cat("\n===============================================\n")
cat("=== PARAMETER RECOVERY QUALITY ===\n")
cat("===============================================\n")
cat("Criterion: r >= 0.7 indicates good recovery\n\n")

# Calculate correlations between true and inferred
recovery_cors <- data.frame(
  Parameter = c("w", "A", "theta", "a"),
  Correlation = c(
    cor(all_true_mu_w, all_infer_mu_w),
    cor(all_true_mu_A, all_infer_mu_A),
    cor(all_true_mu_theta, all_infer_mu_theta),
    cor(all_true_mu_a, all_infer_mu_a)
  ),
  stringsAsFactors = FALSE
)

print(recovery_cors)

# Count parameters with good recovery
n_good <- sum(recovery_cors$Correlation >= 0.7)
cat("\nSummary:", n_good, "out of 4 parameters show good recovery (r >= 0.7)\n")
cat("\nSaved: pvl_recovery_quality_correlations.png\n")

# ===== SIGMA RECOVERY QUALITY =====
cat("\n===============================================\n")
cat("=== SIGMA (BETWEEN-SUBJECT VARIABILITY) RECOVERY ===\n")
cat("===============================================\n")
cat("Converting lambda (precision) to sigma (SD): sigma = 1/sqrt(lambda)\n\n")

# Convert inferred lambda (precision) to sigma (SD)
infer_sigma_w <- 1/sqrt(all_infer_lambda_w)
infer_sigma_A <- 1/sqrt(all_infer_lambda_A)
infer_sigma_theta <- 1/sqrt(all_infer_lambda_theta)
infer_sigma_a <- 1/sqrt(all_infer_lambda_a)

# True values are already stored as sigma (SD)
sigma_recovery_cors <- data.frame(
  Parameter = c("sigma_w", "sigma_A", "sigma_theta", "sigma_a"),
  Correlation = c(
    cor(all_true_lambda_w, infer_sigma_w),
    cor(all_true_lambda_A, infer_sigma_A),
    cor(all_true_lambda_theta, infer_sigma_theta),
    cor(all_true_lambda_a, infer_sigma_a)
  ),
  stringsAsFactors = FALSE
)

print(sigma_recovery_cors)

n_good_sigma <- sum(sigma_recovery_cors$Correlation >= 0.7)
cat("\nSummary:", n_good_sigma, "out of 4 sigma parameters show good recovery (r >= 0.7)\n")

# ===== 1. CORRELATION MATRIX: TRUE PARAMETERS =====
cat("\n===============================================\n")
cat("=== CORRELATIONS BETWEEN TRUE GROUP-LEVEL PARAMETERS ===\n")
cat("===============================================\n")

true_params <- data.frame(
  w = all_true_mu_w,
  A = all_true_mu_A,
  theta = all_true_mu_theta,
  a = all_true_mu_a
)

cor_true <- cor(true_params)
print(round(cor_true, 3))
cat("Saved: pvl_correlation_true_parameters.png\n")

# ===== 2. CORRELATION MATRIX: INFERRED PARAMETERS =====
cat("\n===============================================\n")
cat("=== CORRELATIONS BETWEEN INFERRED GROUP-LEVEL PARAMETERS ===\n")
cat("===============================================\n")

infer_params <- data.frame(
  w = all_infer_mu_w,
  A = all_infer_mu_A,
  theta = all_infer_mu_theta,
  a = all_infer_mu_a
)

cor_infer <- cor(infer_params)
print(round(cor_infer, 3))
cat("Saved: pvl_correlation_inferred_parameters.png\n")

# ===== 3. INVESTIGATING PARAMETER TRADE-OFFS =====
cat("\n===============================================\n")
cat("=== INVESTIGATING PARAMETER TRADE-OFFS ===\n")
cat("===============================================\n")

# Check all pairwise correlations in true parameters
cat("\nTrue parameter correlations:\n")
print(round(cor_true, 3))

# Check all pairwise correlations in inferred parameters
cat("\nInferred parameter correlations:\n")
print(round(cor_infer, 3))

# Visualize strongest relationships (where |r| > 0.3)
cat("\n--- Identifying strong correlations (|r| > 0.3) ---\n")
for (i in 1:3) {
  for (j in (i+1):4) {
    if (abs(cor_true[i,j]) > 0.3) {
      cat(sprintf("True: %s vs %s: r = %.3f\n", 
                  colnames(cor_true)[i], colnames(cor_true)[j], cor_true[i,j]))
    }
    if (abs(cor_infer[i,j]) > 0.3) {
      cat(sprintf("Inferred: %s vs %s: r = %.3f\n", 
                  colnames(cor_infer)[i], colnames(cor_infer)[j], cor_infer[i,j]))
    }
  }
}

# ===== 4. RECOVERY ERRORS: Which parameters cause problems? =====
cat("\n===============================================\n")
cat("=== RECOVERY ERRORS ANALYSIS ===\n")
cat("===============================================\n")

# Calculate recovery errors
errors <- data.frame(
  w_error = all_infer_mu_w - all_true_mu_w,
  A_error = all_infer_mu_A - all_true_mu_A,
  theta_error = all_infer_mu_theta - all_true_mu_theta,
  a_error = all_infer_mu_a - all_true_mu_a
)

# Summary statistics
cat("\nMean Absolute Errors:\n")
mae <- colMeans(abs(errors))
print(round(mae, 4))

cat("\nRoot Mean Squared Errors:\n")
rmse <- sqrt(colMeans(errors^2))
print(round(rmse, 4))

cat("\nMean Bias (positive = overestimation):\n")
bias <- colMeans(errors)
print(round(bias, 4))

# Check if errors correlate with true values (systematic bias)
cat("\n--- BIAS CHECK: Do errors depend on true values? ---\n")
bias_cors <- data.frame(
  Parameter = c("w", "A", "theta", "a"),
  Error_vs_True = c(
    cor(errors$w_error, all_true_mu_w),
    cor(errors$A_error, all_true_mu_A),
    cor(errors$theta_error, all_true_mu_theta),
    cor(errors$a_error, all_true_mu_a)
  )
)
print(bias_cors)
cat("\nInterpretation:\n")
cat("  Negative correlation = underestimation bias at high values\n")
cat("  Positive correlation = overestimation bias at high values\n")
cat("Saved: pvl_recovery_error_patterns.png\n")

# ===== 5. FULL CORRELATION MATRIX: TRUE vs INFERRED =====
cat("\n===============================================\n")
cat("=== FULL CROSS-PARAMETER CORRELATION MATRIX ===\n")
cat("===============================================\n")

all_params <- data.frame(
  true_w = all_true_mu_w,
  true_A = all_true_mu_A,
  true_theta = all_true_mu_theta,
  true_a = all_true_mu_a,
  infer_w = all_infer_mu_w,
  infer_A = all_infer_mu_A,
  infer_theta = all_infer_mu_theta,
  infer_a = all_infer_mu_a
)

cor_full <- cor(all_params)

# Extract cross-correlations (true vs inferred)
cross_cor <- cor_full[1:4, 5:8]
rownames(cross_cor) <- c("true_w", "true_A", "true_theta", "true_a")
colnames(cross_cor) <- c("infer_w", "infer_A", "infer_theta", "infer_a")

cat("\nCross-parameter correlations (rows=true, cols=inferred):\n")
print(round(cross_cor, 3))
cat("Saved: pvl_cross_parameter_correlations.png\n")

# Highlight problematic cross-correlations
cat("\n--- PROBLEMATIC CROSS-CORRELATIONS (|r| > 0.3 off-diagonal) ---\n")
cat("These indicate confusability between parameters:\n\n")
found_issues <- FALSE
for (i in 1:4) {
  for (j in 1:4) {
    if (i != j && abs(cross_cor[i, j]) > 0.3) {
      cat(sprintf("  %s <-> %s: r = %6.3f\n", 
                  rownames(cross_cor)[i], 
                  colnames(cross_cor)[j], 
                  cross_cor[i, j]))
      found_issues <- TRUE
    }
  }
}
if (!found_issues) {
  cat("  None found - parameters are distinguishable!\n")
}

# ===== FINAL SUMMARY =====
cat("\n===============================================\n")
cat("=== SUMMARY ===\n")
cat("===============================================\n")

cat("\nParameter Recovery (r between true and inferred):\n")
print(recovery_cors)

cat("\nSigma Recovery:\n")
print(sigma_recovery_cors)

cat("\n===============================================\n")
cat("=== FILES GENERATED ===\n")
cat("===============================================\n")
cat("  1. pvl_recovery_quality_correlations.png - Main recovery plot\n")
cat("  2. pvl_correlation_true_parameters.png\n")
cat("  3. pvl_correlation_inferred_parameters.png\n")
cat("  4. pvl_recovery_error_patterns.png\n")
cat("  5. pvl_cross_parameter_correlations.png\n")
cat("  6. recovery_analysis_detailed.txt - This file\n")
cat("\n===============================================\n")

# Close sink before creating plots
sink()

# ===== NOW CREATE ALL PLOTS =====
cat("\n===============================================\n")
cat("Creating visualizations...\n")
cat("===============================================\n")

# 1. Recovery quality plot (ggplot2)
p <- ggplot(recovery_cors, aes(x = Parameter, y = Correlation, fill = Correlation >= 0.7)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", size = 1.2) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), 
            vjust = ifelse(recovery_cors$Correlation >= 0, -0.5, 1.5),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
                    labels = c("TRUE" = "Good (>=0.7)", "FALSE" = "Poor (<0.7)"),
                    name = "Recovery Quality") +
  ylim(-1, 1.1) +
  labs(title = "PVL Parameter Recovery Quality (True vs Inferred Correlations)",
       subtitle = "Red dashed line = 0.7 threshold for acceptable recovery",
       y = "Pearson Correlation (r)", x = "Parameter") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 12))

ggsave(file.path(output_dir, "pvl_recovery_quality_correlations.png"), 
       plot = p, width = 8, height = 5, dpi = 100)
cat("✓ Saved: pvl_recovery_quality_correlations.png\n")

# 2. True parameters correlation matrix
png(file.path(output_dir, "pvl_correlation_true_parameters.png"), width = 800, height = 800)
corrplot(cor_true, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Correlations: True Group-Level Parameters (PVL)",
         mar = c(0,0,2,0))
dev.off()
cat("✓ Saved: pvl_correlation_true_parameters.png\n")

# 3. Inferred parameters correlation matrix
png(file.path(output_dir, "pvl_correlation_inferred_parameters.png"), width = 800, height = 800)
corrplot(cor_infer, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Correlations: Inferred Group-Level Parameters (PVL)",
         mar = c(0,0,2,0))
dev.off()
cat("✓ Saved: pvl_correlation_inferred_parameters.png\n")

# 4. Error patterns
png(file.path(output_dir, "pvl_recovery_error_patterns.png"), width = 1200, height = 600)
par(mfrow = c(2, 2))

param_names <- c("w", "A", "theta", "a")
for (param in param_names) {
  true_vals <- get(paste0("all_true_mu_", param))
  errors_vals <- errors[[paste0(param, "_error")]]
  
  plot(true_vals, errors_vals,
       xlab = paste("True mu_", param, sep=""),
       ylab = "Recovery Error (Inferred - True)",
       main = paste(param, "\nBias r =", round(cor(errors_vals, true_vals), 3)),
       pch = 19, col = "blue", cex = 1.2)
  abline(h = 0, col = "red", lty = 2, lwd = 2)
  abline(lm(errors_vals ~ true_vals), col = "darkgreen", lwd = 2)
}

dev.off()
cat("✓ Saved: pvl_recovery_error_patterns.png\n")

# 5. Cross-parameter correlations
png(file.path(output_dir, "pvl_cross_parameter_correlations.png"), width = 900, height = 700)
corrplot(cross_cor, method = "color", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Cross-Parameter Correlations (PVL)\n(True vs Inferred)",
         mar = c(0,0,2,0),
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()
cat("✓ Saved: pvl_cross_parameter_correlations.png\n")

# Print final summary to console
cat("\n===============================================\n")
cat("=== PARAMETER RECOVERY ANALYSIS COMPLETE ===\n")
cat("===============================================\n")
cat("\nParameter Recovery (r between true and inferred):\n")
print(recovery_cors)
cat("\nSigma Recovery:\n")
print(sigma_recovery_cors)
cat("\nAll detailed results saved to:", analysis_file, "\n")
cat("All plots saved to:", output_dir, "\n")
cat("===============================================\n")