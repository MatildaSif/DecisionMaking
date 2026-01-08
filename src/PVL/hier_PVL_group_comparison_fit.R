# PVL Group Comparison Analysis
library(pacman)
pacman::p_load(extraDistr, R2jags, parallel, ggpubr, rstatix, gridExtra, 
               polspline, truncnorm, ggplot2, tidyr, dplyr)

set.seed(1702)

setwd('/work/JoMat/DecisionMaking/')

# Output directory
output_dir <- "src/PVL/outputs/group_comparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define MPD function
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------Load data----
ctr_data <- read.table("data/IGTdata_healthy_control_clean.txt", header=TRUE)
her_data <- read.table("data/IGTdata_heroin_clean.txt", header=TRUE)
amp_data <- read.table("data/IGTdata_amphetamine_clean.txt", header=TRUE)

#----------Prepare data for JAGS----
# Identify and count unique subject IDs
subIDs_ctr <- unique(ctr_data$subjID)
nsubs_ctr <- length(subIDs_ctr)

subIDs_her <- unique(her_data$subjID)
nsubs_her <- length(subIDs_her)

subIDs_amp <- unique(amp_data$subjID)
nsubs_amp <- length(subIDs_amp)

ntrials_max <- 100

# Get raw data
x_raw_ctr <- ctr_data$deck
X_raw_ctr <- (ctr_data$gain + ctr_data$loss)

x_raw_her <- her_data$deck
X_raw_her <- (her_data$gain + her_data$loss)

x_raw_amp <- amp_data$deck
X_raw_amp <- (amp_data$gain + amp_data$loss)

#--- Assign choices and outcomes in trial x sub matrix ---

# Empty arrays to fill
ntrials_ctr <- array(0, c(nsubs_ctr))
x_ctr <- array(0, c(nsubs_ctr, ntrials_max))
X_ctr <- array(0, c(nsubs_ctr, ntrials_max))

ntrials_her <- array(0, c(nsubs_her))
x_her <- array(0, c(nsubs_her, ntrials_max))
X_her <- array(0, c(nsubs_her, ntrials_max))

ntrials_amp <- array(0, c(nsubs_amp))
x_amp <- array(0, c(nsubs_amp, ntrials_max))
X_amp <- array(0, c(nsubs_amp, ntrials_max))

# Make control data matrices
for (s in 1:nsubs_ctr) {
  ntrials_ctr[s] <- length(x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]])
  
  x_sub_ctr <- x_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(x_sub_ctr) <- ntrials_max
  
  X_sub_ctr <- X_raw_ctr[ctr_data$subjID == subIDs_ctr[s]] 
  length(X_sub_ctr) <- ntrials_max
  
  x_ctr[s,] <- x_sub_ctr
  X_ctr[s,] <- X_sub_ctr
}

# Make heroin data matrices
for (s in 1:nsubs_her) {
  ntrials_her[s] <- length(x_raw_her[her_data$subjID == subIDs_her[s]])
  
  x_sub_her <- x_raw_her[her_data$subjID == subIDs_her[s]] 
  length(x_sub_her) <- ntrials_max
  
  X_sub_her <- X_raw_her[her_data$subjID == subIDs_her[s]] 
  length(X_sub_her) <- ntrials_max
  
  x_her[s,] <- x_sub_her
  X_her[s,] <- X_sub_her
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

#----------Run PVL model comparison----
cat("\n=== Running PVL Group Comparison Model ===\n")
cat("Controls:", nsubs_ctr, "subjects\n")
cat("Heroin:", nsubs_her, "subjects\n")
cat("Amphetamine:", nsubs_amp, "subjects\n\n")

# Set up JAGS data
data <- list("x_ctr", "X_ctr", "ntrials_ctr", "nsubs_ctr",
             "x_her", "X_her", "ntrials_her", "nsubs_her",
             "x_amp", "X_amp", "ntrials_amp", "nsubs_amp")

# Parameters to track
params <- c("mu_w_ctr", "mu_A_ctr", "mu_theta_ctr", "mu_a_ctr",
            "mu_w_her", "mu_A_her", "mu_theta_her", "mu_a_her",
            "mu_w_amp", "mu_A_amp", "mu_theta_amp", "mu_a_amp",
            "diff_w_ctr_her", "diff_A_ctr_her", "diff_theta_ctr_her", "diff_a_ctr_her",
            "diff_w_ctr_amp", "diff_A_ctr_amp", "diff_theta_ctr_amp", "diff_a_ctr_amp",
            "diff_w_her_amp", "diff_A_her_amp", "diff_theta_her_amp", "diff_a_her_amp",
            "p_ctr", "p_her", "p_amp")

# Run JAGS
start_time <- Sys.time()

samples <- jags.parallel(data, inits=NULL, params,
                         model.file="src/PVL/hier_PVL_group.txt",
                         n.chains=3, 
                         n.iter=3000, 
                         n.burnin=1000, 
                         n.thin=1,
                         n.cluster=4)

end_time <- Sys.time()

cat("\nModel fitting completed in:", 
    round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")


#----------Save results----
save(samples, file = file.path(output_dir, "PVL_group_comparison_results.RData"))
cat("\nSaved:", file.path(output_dir, "PVL_group_comparison_results.RData"), "\n")


#----------Extract and summarize results----
Y <- samples$BUGSoutput$sims.list

# Group-level parameters
cat("\n=== GROUP-LEVEL PARAMETER ESTIMATES ===\n")

cat("\n--- CONTROLS ---\n")
cat("w (loss aversion):       ", round(MPD(Y$mu_w_ctr), 3), "\n")
cat("A (outcome sensitivity): ", round(MPD(Y$mu_A_ctr), 3), "\n")
cat("theta (consistency):     ", round(MPD(Y$mu_theta_ctr), 3), "\n")
cat("a (learning rate):       ", round(MPD(Y$mu_a_ctr), 3), "\n")

cat("\n--- HEROIN USERS ---\n")
cat("w (loss aversion):       ", round(MPD(Y$mu_w_her), 3), "\n")
cat("A (outcome sensitivity): ", round(MPD(Y$mu_A_her), 3), "\n")
cat("theta (consistency):     ", round(MPD(Y$mu_theta_her), 3), "\n")
cat("a (learning rate):       ", round(MPD(Y$mu_a_her), 3), "\n")

cat("\n--- AMPHETAMINE USERS ---\n")
cat("w (loss aversion):       ", round(MPD(Y$mu_w_amp), 3), "\n")
cat("A (outcome sensitivity): ", round(MPD(Y$mu_A_amp), 3), "\n")
cat("theta (consistency):     ", round(MPD(Y$mu_theta_amp), 3), "\n")
cat("a (learning rate):       ", round(MPD(Y$mu_a_amp), 3), "\n")

#----------Visualize group differences----

# Function to create comparison plot
compare_plot <- function(param_name, ctr, her, amp, xlim_range) {
  plot(density(ctr), xlim=xlim_range, main=param_name, 
       col="lightgreen", lwd=2, xlab="Parameter value", ylab="Density")
  lines(density(her), col="lightblue", lwd=2)
  lines(density(amp), col="salmon", lwd=2)
  legend("topright", legend=c("Controls", "Heroin", "Amphetamine"), 
         col=c("lightgreen", "lightblue", "salmon"), lwd=2)
  abline(v=0, lty=2, col="black")
}

# Create comparison plots
png(file.path(output_dir, "PVL_group_comparison_parameters.png"), 
    width=1200, height=800)
par(mfrow=c(2,2))

compare_plot("w (loss aversion)", 
             Y$mu_w_ctr, Y$mu_w_her, Y$mu_w_amp, c(0, 3))
compare_plot("A (outcome sensitivity)", 
             Y$mu_A_ctr, Y$mu_A_her, Y$mu_A_amp, c(-1, 2))
compare_plot("theta (consistency)", 
             Y$mu_theta_ctr, Y$mu_theta_her, Y$mu_theta_amp, c(0, 5))
compare_plot("a (learning rate)", 
             Y$mu_a_ctr, Y$mu_a_her, Y$mu_a_amp, c(0, 1))

dev.off()
cat("\nSaved:", file.path(output_dir, "PVL_group_comparison_parameters.png"), "\n")

#----------Prior vs Posterior Plots----

plot_prior_posterior <- function(posterior_samples, prior_mean, prior_sd, 
                                 lower_bound = -Inf, upper_bound = Inf,
                                 param_name, xlim_range) {
  
  x_seq <- seq(xlim_range[1], xlim_range[2], length.out = 1000)
  
  if (is.finite(lower_bound) && is.finite(upper_bound)) {
    prior_density <- dtruncnorm(x_seq, a = lower_bound, b = upper_bound, 
                                mean = prior_mean, sd = prior_sd)
  } else if (is.finite(lower_bound)) {
    prior_density <- dtruncnorm(x_seq, a = lower_bound, b = Inf, 
                                mean = prior_mean, sd = prior_sd)
  } else {
    prior_density <- dnorm(x_seq, mean = prior_mean, sd = prior_sd)
  }
  
  plot(density(posterior_samples), main = param_name, 
       xlab = "Parameter value", ylab = "Density",
       col = "red", lwd = 2, xlim = xlim_range, ylim = c(0, max(density(posterior_samples)$y) * 1.1))
  polygon(density(posterior_samples), col = rgb(0, 0, 1, 0.3), border = NA)
  
  lines(x_seq, prior_density, col = "black", lwd = 2, lty = 2)
  
  legend("topright", legend = c("Posterior", "Prior"), 
         col = c("red", "black"), lwd = 2, lty = c(1, 2), bty = "n")
}

# Prior-posterior plots for CONTROLS
png(file.path(output_dir, "PVL_prior_posterior_controls.png"), 
    width = 1200, height = 800)
par(mfrow = c(2, 2))

plot_prior_posterior(Y$mu_w_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = Inf,
                     param_name = "w (Controls)", xlim_range = c(0, 3))

plot_prior_posterior(Y$mu_A_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = -Inf, upper_bound = Inf,
                     param_name = "A (Controls)", xlim_range = c(-2, 2))

plot_prior_posterior(Y$mu_theta_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = Inf,
                     param_name = "theta (Controls)", xlim_range = c(0, 5))

plot_prior_posterior(Y$mu_a_ctr, prior_mean = 0, prior_sd = 1, 
                     lower_bound = 0, upper_bound = 1,
                     param_name = "a (Controls)", xlim_range = c(0, 1))

dev.off()
cat("Saved:", file.path(output_dir, "PVL_prior_posterior_controls.png"), "\n")

# Repeat for HEROIN
png(file.path(output_dir, "PVL_prior_posterior_heroin.png"), 
    width = 1200, height = 800)
par(mfrow = c(2, 2))

plot_prior_posterior(Y$mu_w_her, 0, 1, 0, Inf, "w (Heroin)", c(0, 3))
plot_prior_posterior(Y$mu_A_her, 0, 1, -Inf, Inf, "A (Heroin)", c(-2, 2))
plot_prior_posterior(Y$mu_theta_her, 0, 1, 0, Inf, "theta (Heroin)", c(0, 5))
plot_prior_posterior(Y$mu_a_her, 0, 1, 0, 1, "a (Heroin)", c(0, 1))

dev.off()
cat("Saved:", file.path(output_dir, "PVL_prior_posterior_heroin.png"), "\n")

# Repeat for AMPHETAMINE
png(file.path(output_dir, "PVL_prior_posterior_amphetamine.png"), 
    width = 1200, height = 800)
par(mfrow = c(2, 2))

plot_prior_posterior(Y$mu_w_amp, 0, 1, 0, Inf, "w (Amphetamine)", c(0, 3))
plot_prior_posterior(Y$mu_A_amp, 0, 1, -Inf, Inf, "A (Amphetamine)", c(-2, 2))
plot_prior_posterior(Y$mu_theta_amp, 0, 1, 0, Inf, "theta (Amphetamine)", c(0, 5))
plot_prior_posterior(Y$mu_a_amp, 0, 1, 0, 1, "a (Amphetamine)", c(0, 1))

dev.off()
cat("Saved:", file.path(output_dir, "PVL_prior_posterior_amphetamine.png"), "\n")

#----------Statistical Testing----

# Calculate net outcomes
calculate_net_outcome <- function(X_mat, ntrials_vec, nsubs) {
  net_outcomes <- rep(NA, nsubs)
  for (s in 1:nsubs) {
    outcomes <- X_mat[s, 1:ntrials_vec[s]]
    outcomes <- na.omit(outcomes)
    net_outcomes[s] <- sum(outcomes)
  }
  return(net_outcomes)
}

net_score_ctr <- calculate_net_outcome(X_ctr, ntrials_ctr, nsubs_ctr)
net_score_her <- calculate_net_outcome(X_her, ntrials_her, nsubs_her)
net_score_amp <- calculate_net_outcome(X_amp, ntrials_amp, nsubs_amp)

# Create dataframe
df_net_score <- data.frame(
  value = c(net_score_ctr, net_score_her, net_score_amp),
  group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                     c(nsubs_ctr, nsubs_her, nsubs_amp)),
                 levels = c("Heroin", "Controls", "Amphetamine"))
)

# Function for complete statistical analysis
run_stats <- function(df, measure_name) {
  cat("\n===============================================\n")
  cat(measure_name, "\n")
  cat("===============================================\n")
  
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
  
  anova_result <- aov(value ~ group, data = df)
  anova_summary <- summary(anova_result)
  
  cat("\nOne-way ANOVA:\n")
  print(anova_summary)
  
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]
  
  cat("\nF(", anova_summary[[1]]$Df[1], ", ", anova_summary[[1]]$Df[2], 
      ") = ", round(f_stat, 3), ", p = ", round(p_value, 4), "\n", sep = "")
  
  if (p_value < 0.05) {
    cat("\n*** ANOVA significant - Running post-hoc tests ***\n")
    posthoc <- TukeyHSD(anova_result)
    cat("\nTukey's HSD Post-hoc Comparisons:\n")
    print(posthoc)
    
    posthoc_df <- as.data.frame(posthoc$group)
    posthoc_df$comparison <- rownames(posthoc_df)
    posthoc_df$significant <- posthoc_df$`p adj` < 0.05
    
    return(list(anova_p = p_value, posthoc = posthoc_df, significant = TRUE))
  } else {
    cat("\n*** ANOVA not significant - No post-hoc tests needed ***\n")
    return(list(anova_p = p_value, posthoc = NULL, significant = FALSE))
  }
}

# Run statistics
stats_net <- run_stats(df_net_score, "NET SCORE")

# Plot with significance
plot_with_significance <- function(df, stats_result, y_label, title, y_limits = NULL) {
  
  summary_stats <- df %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n()),
      .groups = 'drop'
    )
  
  colors <- c("Heroin" = "lightblue", "Controls" = "lightgreen", "Amphetamine" = "salmon")
  
  p <- ggplot(df, aes(x = group, y = value, color = group, fill = group)) +
    geom_bar(data = summary_stats, aes(x = group, y = mean), 
             stat = "identity", alpha = 0.6, width = 0.7) +
    geom_errorbar(data = summary_stats, 
                  aes(x = group, y = mean, ymin = mean - sd, ymax = mean + sd),
                  width = 0.2, color = "black", size = 0.8) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.6) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(title = title, x = "", y = y_label) +
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
  
  if (stats_result$significant && !is.null(stats_result$posthoc)) {
    y_max <- max(df$value, na.rm = TRUE)
    y_min <- min(df$value, na.rm = TRUE)
    y_range <- y_max - y_min
    
    sig_comparisons <- stats_result$posthoc[stats_result$posthoc$significant, ]
    
    if (nrow(sig_comparisons) > 0) {
      y_start <- y_max + y_range * 0.05
      
      for (i in 1:nrow(sig_comparisons)) {
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

# Create plot
p_net <- plot_with_significance(df_net_score, stats_net, "Net Score", "Net Score")

ggsave(file.path(output_dir, "PVL_Net_Score_Comparison.png"), 
       p_net, width = 6, height = 5, dpi = 300)

cat("\n\n=== PLOT SAVED ===\n")
cat("File:", file.path(output_dir, "PVL_Net_Score_Comparison.png"), "\n")


cat("\n✅ PVL group comparison analysis complete.\n")