# Load the recovery results
load("VSE_recovery_results_aspaper.RData")

# Install/load required packages
if (!require("corrplot")) install.packages("corrplot")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
library(corrplot)
library(ggplot2)
library(reshape2)

# ===== PARAMETER RECOVERY QUALITY (PRIMARY ANALYSIS) =====
cat("\n===============================================\n")
cat("=== PARAMETER RECOVERY QUALITY ===\n")
cat("===============================================\n")
cat("Criterion: r >= 0.7 indicates good recovery\n\n")

# Calculate correlations between true and inferred
recovery_cors <- data.frame(
  Parameter = c("theta", "delta", "alpha", "phi", "c"),
  Correlation = c(
    cor(true_mu_theta, infer_mu_theta),
    cor(true_mu_delta, infer_mu_delta),
    cor(true_mu_alpha, infer_mu_alpha),
    cor(true_mu_phi, infer_mu_phi),
    cor(true_mu_c, infer_mu_c)
  ),
  stringsAsFactors = FALSE
)

print(recovery_cors)

# Count parameters with good recovery
n_good <- sum(recovery_cors$Correlation >= 0.7)
cat("\nSummary:", n_good, "out of 5 parameters show good recovery (r >= 0.7)\n")

# Visualize recovery quality
png("recovery_quality_correlations.png", width = 800, height = 500)
ggplot(recovery_cors, aes(x = Parameter, y = Correlation, fill = Correlation >= 0.7)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", size = 1.2) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), 
            vjust = ifelse(recovery_cors$Correlation >= 0, -0.5, 1.5),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
                    labels = c("TRUE" = "Good (>=0.7)", "FALSE" = "Poor (<0.7)"),
                    name = "Recovery Quality") +
  ylim(-1, 1.1) +
  labs(title = "Parameter Recovery Quality (True vs Inferred Correlations)",
       subtitle = "Red dashed line = 0.7 threshold for acceptable recovery",
       y = "Pearson Correlation (r)", x = "Parameter") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 12))
dev.off()
cat("\nSaved: recovery_quality_correlations.png\n")

# ===== SIGMA RECOVERY QUALITY =====
cat("\n===============================================\n")
cat("=== SIGMA (BETWEEN-SUBJECT VARIABILITY) RECOVERY ===\n")
cat("===============================================\n")

sigma_recovery_cors <- data.frame(
  Parameter = c("sigma_theta", "sigma_delta", "sigma_alpha", "sigma_phi", "sigma_c"),
  Correlation = c(
    cor(true_sigma_theta, infer_sigma_theta),
    cor(true_sigma_delta, infer_sigma_delta),
    cor(true_sigma_alpha, infer_sigma_alpha),
    cor(true_sigma_phi, infer_sigma_phi),
    cor(true_sigma_c, infer_sigma_c)
  ),
  stringsAsFactors = FALSE
)

print(sigma_recovery_cors)

n_good_sigma <- sum(sigma_recovery_cors$Correlation >= 0.7)
cat("\nSummary:", n_good_sigma, "out of 5 sigma parameters show good recovery (r >= 0.7)\n")

# ===== 1. CORRELATION MATRIX: TRUE PARAMETERS =====
cat("\n===============================================\n")
cat("=== CORRELATIONS BETWEEN TRUE GROUP-LEVEL PARAMETERS ===\n")
cat("===============================================\n")

true_params <- data.frame(
  theta = true_mu_theta,
  delta = true_mu_delta,
  alpha = true_mu_alpha,
  phi = true_mu_phi,
  c = true_mu_c
)

cor_true <- cor(true_params)
print(round(cor_true, 3))

# Visualize
png("correlation_true_parameters.png", width = 800, height = 800)
corrplot(cor_true, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Correlations: True Group-Level Parameters",
         mar = c(0,0,2,0))
dev.off()
cat("Saved: correlation_true_parameters.png\n")

# ===== 2. CORRELATION MATRIX: INFERRED PARAMETERS =====
cat("\n===============================================\n")
cat("=== CORRELATIONS BETWEEN INFERRED GROUP-LEVEL PARAMETERS ===\n")
cat("===============================================\n")

infer_params <- data.frame(
  theta = infer_mu_theta,
  delta = infer_mu_delta,
  alpha = infer_mu_alpha,
  phi = infer_mu_phi,
  c = infer_mu_c
)

cor_infer <- cor(infer_params)
print(round(cor_infer, 3))

# Visualize
png("correlation_inferred_parameters.png", width = 800, height = 800)
corrplot(cor_infer, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Correlations: Inferred Group-Level Parameters",
         mar = c(0,0,2,0))
dev.off()
cat("Saved: correlation_inferred_parameters.png\n")

# ===== 3. INVESTIGATING THETA-DELTA TRADE-OFF =====
cat("\n===============================================\n")
cat("=== INVESTIGATING THETA-DELTA TRADE-OFF ===\n")
cat("===============================================\n")

# Theta-Delta relationship
cor_theta_delta_true <- cor(true_mu_theta, true_mu_delta)
cat("True theta-delta correlation:", round(cor_theta_delta_true, 3), "\n")

cor_theta_delta_infer <- cor(infer_mu_theta, infer_mu_delta)
cat("Inferred theta-delta correlation:", round(cor_theta_delta_infer, 3), "\n")

# Check cross-correlation
cor_true_theta_infer_delta <- cor(true_mu_theta, infer_mu_delta)
cor_true_delta_infer_theta <- cor(true_mu_delta, infer_mu_theta)
cat("True theta vs Inferred delta:", round(cor_true_theta_infer_delta, 3), "\n")
cat("True delta vs Inferred theta:", round(cor_true_delta_infer_theta, 3), "\n")

# Visualize theta-delta relationship
png("theta_delta_relationship.png", width = 1000, height = 400)
par(mfrow = c(1, 2))

plot(true_mu_theta, true_mu_delta, 
     xlab = "True mu_theta", ylab = "True mu_delta",
     main = paste("True Parameters\nr =", round(cor_theta_delta_true, 3)),
     pch = 19, col = "blue", cex = 1.2)
abline(lm(true_mu_delta ~ true_mu_theta), col = "red", lwd = 2)

plot(infer_mu_theta, infer_mu_delta, 
     xlab = "Inferred mu_theta", ylab = "Inferred mu_delta",
     main = paste("Inferred Parameters\nr =", round(cor_theta_delta_infer, 3)),
     pch = 19, col = "darkgreen", cex = 1.2)
abline(lm(infer_mu_delta ~ infer_mu_theta), col = "red", lwd = 2)

dev.off()
cat("Saved: theta_delta_relationship.png\n")

# ===== 4. RECOVERY ERRORS: Which parameters cause problems? =====
cat("\n===============================================\n")
cat("=== RECOVERY ERRORS ANALYSIS ===\n")
cat("===============================================\n")

# Calculate recovery errors
errors <- data.frame(
  theta_error = infer_mu_theta - true_mu_theta,
  delta_error = infer_mu_delta - true_mu_delta,
  alpha_error = infer_mu_alpha - true_mu_alpha,
  phi_error = infer_mu_phi - true_mu_phi,
  c_error = infer_mu_c - true_mu_c
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
  Parameter = c("theta", "delta", "alpha", "phi", "c"),
  Error_vs_True = c(
    cor(errors$theta_error, true_mu_theta),
    cor(errors$delta_error, true_mu_delta),
    cor(errors$alpha_error, true_mu_alpha),
    cor(errors$phi_error, true_mu_phi),
    cor(errors$c_error, true_mu_c)
  )
)
print(bias_cors)
cat("\nInterpretation:\n")
cat("  Negative correlation = underestimation bias at high values\n")
cat("  Positive correlation = overestimation bias at high values\n")

# Visualize error patterns
png("recovery_error_patterns.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

for (param in c("theta", "delta", "alpha", "phi", "c")) {
  true_vals <- get(paste0("true_mu_", param))
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
cat("Saved: recovery_error_patterns.png\n")

# ===== 5. FULL CORRELATION MATRIX: TRUE vs INFERRED =====
cat("\n===============================================\n")
cat("=== FULL CROSS-PARAMETER CORRELATION MATRIX ===\n")
cat("===============================================\n")

all_params <- data.frame(
  true_theta = true_mu_theta,
  true_delta = true_mu_delta,
  true_alpha = true_mu_alpha,
  true_phi = true_mu_phi,
  true_c = true_mu_c,
  infer_theta = infer_mu_theta,
  infer_delta = infer_mu_delta,
  infer_alpha = infer_mu_alpha,
  infer_phi = infer_mu_phi,
  infer_c = infer_mu_c
)

cor_full <- cor(all_params)

# Extract cross-correlations (true vs inferred)
cross_cor <- cor_full[1:5, 6:10]
rownames(cross_cor) <- c("true_theta", "true_delta", "true_alpha", "true_phi", "true_c")
colnames(cross_cor) <- c("infer_theta", "infer_delta", "infer_alpha", "infer_phi", "infer_c")

cat("\nCross-parameter correlations (rows=true, cols=inferred):\n")
print(round(cross_cor, 3))

png("cross_parameter_correlations.png", width = 900, height = 700)
corrplot(cross_cor, method = "color", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Cross-Parameter Correlations\n(True vs Inferred)",
         mar = c(0,0,2,0),
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()
cat("Saved: cross_parameter_correlations.png\n")

# Highlight problematic cross-correlations
cat("\n--- PROBLEMATIC CROSS-CORRELATIONS (|r| > 0.3 off-diagonal) ---\n")
cat("These indicate confusability between parameters:\n\n")
found_issues <- FALSE
for (i in 1:5) {
  for (j in 1:5) {
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
cat("  1. recovery_quality_correlations.png - Main recovery plot\n")
cat("  2. correlation_true_parameters.png\n")
cat("  3. correlation_inferred_parameters.png\n")
cat("  4. theta_delta_relationship.png\n")
cat("  5. recovery_error_patterns.png\n")
cat("  6. cross_parameter_correlations.png\n")
cat("\n===============================================\n")