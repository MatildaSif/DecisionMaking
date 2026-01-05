library(R2jags)
library(truncnorm)
library(extraDistr)
library(ggplot2)

setwd('/work/JoMat/DecisionMaking/')

# Create output folder structure
output_dir <- "PVL/outputs/parameter_estimation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

source('src/PVL/hier_PVL_sim.R')

# -----------------------------
# Define group
# -----------------------------
GROUP <- "heroin"  # Change to "healthy", "heroin", or "amphetamine"

rm(list = c("samples", "x_all", "ntrials_all", "Y", "subIDs", "nsubs"))

if (GROUP == "healthy") {
  load(file.path(output_dir, "PVL_healthy_posteriors.RData"))
  group_title <- "Healthy Controls"
} else if (GROUP == "heroin") {
  load(file.path(output_dir, "PVL_heroin_posteriors.RData"))
  group_title <- "Heroin Users"
} else if (GROUP == "amphetamine") {
  load(file.path(output_dir, "PVL_amphetamine_posteriors.RData"))
  group_title <- "Amphetamine Users"
}

# -----------------------------
# Recreate payoff structure
# -----------------------------
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

A_R <- rep(bad_r, nstruct)
A_L <- c(rep(bad_freq_l, nstruct*freq), rep(0, nstruct*(1-freq)))
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq)))
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq), rep(0, nstruct*(1-freq)))
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq)))

A <- array(NA, ntrials)
B <- array(NA, ntrials)
C <- array(NA, ntrials)
D <- array(NA, ntrials)

set.seed(1983)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- A_R + sample(A_L)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- B_R + sample(B_L)
  C[(1+(i-1)*nstruct):(i*nstruct)] <- C_R + sample(C_L)
  D[(1+(i-1)*nstruct):(i*nstruct)] <- D_R + sample(D_L)
}

payoff <- cbind(A, B, C, D) / 100

# -----------------------------
# Extract posterior samples
# -----------------------------
Y <- samples$BUGSoutput$sims.list

# Check dimensions
cat("\n=== Checking loaded parameters ===\n")
cat("Group:", group_title, "\n")
cat("Number of subjects:", nsubs, "\n")
cat("Subject-level w dimensions:", dim(Y$w), "\n")
cat("Subject-level A dimensions:", dim(Y$A), "\n")

# Helper function: get MAP
MPD <- function(x) {
  density(x)$x[which.max(density(x)$y)]
}

# -----------------------------
# Subject-level posterior predictive check
# -----------------------------
cat("\n=== Running subject-level PPC ===\n")
n_ppc_samples <- 100  # Number of posterior samples to use per subject

subject_match <- numeric(nsubs)

for (s in 1:nsubs) {
  
  cat("Subject", s, "of", nsubs, "\n")
  
  # Get MAP estimates for THIS SPECIFIC SUBJECT
  w_s <- MPD(Y$w[, s])
  A_s <- MPD(Y$A[, s])
  theta_s <- MPD(Y$theta[, s])
  a_s <- MPD(Y$a[, s])
  
  # Simulate choices using THIS subject's parameters
  # Note: sigmas set to 0 because we're using point estimates
  PVL_sims <- hier_PVL_sim(payoff, 1, ntrials_all[s], 
                           w_s, A_s, a_s, theta_s,
                           0, 0, 0, 0)  # No between-subject variability
  
  predicted_choices <- PVL_sims$x[1, 1:ntrials_all[s]]
  actual_choices <- x_all[s, 1:ntrials_all[s]]
  
  # Calculate match proportion
  subject_match[s] <- sum(predicted_choices == actual_choices, na.rm = TRUE) / ntrials_all[s]
}

# -----------------------------
# Alternative: Sample from posterior uncertainty
# -----------------------------
cat("\n=== Running PPC with posterior uncertainty ===\n")
subject_match_uncertain <- numeric(nsubs)

for (s in 1:nsubs) {
  
  # Sample multiple posterior parameter sets
  n_samples <- min(n_ppc_samples, nrow(Y$w))
  sample_idx <- sample(1:nrow(Y$w), n_samples)
  
  match_samples <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    # Get one posterior sample for this subject
    w_s <- Y$w[sample_idx[i], s]
    A_s <- Y$A[sample_idx[i], s]
    theta_s <- Y$theta[sample_idx[i], s]
    a_s <- Y$a[sample_idx[i], s]
    
    # Simulate
    PVL_sims <- hier_PVL_sim(payoff, 1, ntrials_all[s], 
                             w_s, A_s, a_s, theta_s,
                             0, 0, 0, 0)
    
    predicted_choices <- PVL_sims$x[1, 1:ntrials_all[s]]
    actual_choices <- x_all[s, 1:ntrials_all[s]]
    
    match_samples[i] <- sum(predicted_choices == actual_choices, na.rm = TRUE) / ntrials_all[s]
  }
  
  # Average across posterior samples
  subject_match_uncertain[s] <- mean(match_samples)
}

# Save summary to text file
summary_file <- file.path(output_dir, paste0("PPC_summary_", GROUP, ".txt"))
sink(summary_file)

cat("\n=== SUBJECT-LEVEL POSTERIOR PREDICTIVE CHECK ===\n")
cat("Group:", group_title, "\n\n")

cat("Using MAP estimates:\n")
cat("  Mean accuracy:", round(mean(subject_match), 3), "\n")
cat("  SD:", round(sd(subject_match), 3), "\n")
cat("  Range:", round(min(subject_match), 3), "-", round(max(subject_match), 3), "\n")

cat("\nUsing posterior uncertainty:\n")
cat("  Mean accuracy:", round(mean(subject_match_uncertain), 3), "\n")
cat("  SD:", round(sd(subject_match_uncertain), 3), "\n")
cat("  Range:", round(min(subject_match_uncertain), 3), "-", round(max(subject_match_uncertain), 3), "\n")

cat("\nChance level: 0.25 (25%)\n")

sink()

# Print to console as well
cat("\n=== SUBJECT-LEVEL POSTERIOR PREDICTIVE CHECK ===\n")
cat("Group:", group_title, "\n\n")

cat("Using MAP estimates:\n")
cat("  Mean accuracy:", round(mean(subject_match), 3), "\n")
cat("  SD:", round(sd(subject_match), 3), "\n")
cat("  Range:", round(min(subject_match), 3), "-", round(max(subject_match), 3), "\n")

cat("\nUsing posterior uncertainty:\n")
cat("  Mean accuracy:", round(mean(subject_match_uncertain), 3), "\n")
cat("  SD:", round(sd(subject_match_uncertain), 3), "\n")
cat("  Range:", round(min(subject_match_uncertain), 3), "-", round(max(subject_match_uncertain), 3), "\n")

cat("\nChance level: 0.25 (25%)\n")

# -----------------------------
# Visualization
# -----------------------------
ppc_df <- data.frame(
  subject = factor(1:nsubs),
  MAP = subject_match,
  Posterior = subject_match_uncertain
)

# Reshape for plotting
ppc_long <- tidyr::pivot_longer(ppc_df, cols = c(MAP, Posterior),
                                names_to = "Method", values_to = "Accuracy")

p <- ggplot(ppc_long, aes(x = subject, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed", size = 1) +
  labs(title = paste("Subject-level Posterior Predictive Accuracy:", group_title),
       subtitle = "Red line = chance (25%)",
       x = "Participant",
       y = "Proportion Correctly Predicted Choices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("MAP" = "skyblue", "Posterior" = "lightcoral"))

print(p)

# Save figure
ggsave(filename = file.path(output_dir, paste0("PVL_PPC_subject_", GROUP, ".png")), 
       plot = p, width = 14, height = 6, dpi = 300)
cat("Saved:", paste0("PVL_PPC_subject_", GROUP, ".png\n"))

# Save results
save(subject_match, subject_match_uncertain, ppc_df,
     file = file.path(output_dir, paste0("PVL_PPC_results_", GROUP, ".RData")))
cat("Saved:", paste0("PVL_PPC_results_", GROUP, ".RData\n"))
