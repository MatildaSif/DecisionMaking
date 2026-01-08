# Load required packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

# Set working directory
setwd('/work/JoMat/DecisionMaking/src/VSE')

# Load the saved JAGS output
load("VSE_group_comparison_results.RData")  # Adjust filename as needed
# Assumes the object is named 'samples' or 'Y'

# If samples object exists, extract sims.list
if(exists("samples")) {
  Y <- samples$BUGSoutput$sims.list
}

# Define colors matching your original script
colors <- c("Controls" = "#90EE90",      # lightgreen
            "Heroin" = "#ADD8E6",         # lightblue  
            "Amphetamine" = "#FA8072")    # salmon

# Function to create density plot for each parameter
plot_parameter <- function(ctr_samples, opi_samples, amp_samples, 
                           param_label, xlim_range) {
  
  # Prepare data
  df <- data.frame(
    value = c(ctr_samples, opi_samples, amp_samples),
    group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                       times = c(length(ctr_samples), 
                                 length(opi_samples), 
                                 length(amp_samples))),
                   levels = c("Controls", "Heroin", "Amphetamine"))
  )
  
  # Create plot
  p <- ggplot(df, aes(x = value, fill = group, color = group)) +
    geom_density(alpha = 0.5, size = 1) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(title = param_label,
         x = "",
         y = "Density") +
    xlim(xlim_range) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
  
  return(p)
}

# Create individual plots for each parameter
p_theta <- plot_parameter(Y$mu_theta_ctr, Y$mu_theta_opi, Y$mu_theta_amp,
                          expression(mu[theta]), c(0, 1))

p_delta <- plot_parameter(Y$mu_delta_ctr, Y$mu_delta_opi, Y$mu_delta_amp,
                          expression(mu[delta]), c(0, 1))

p_alpha <- plot_parameter(Y$mu_alpha_ctr, Y$mu_alpha_opi, Y$mu_alpha_amp,
                          expression(mu[alpha]), c(0, 1))

p_phi <- plot_parameter(Y$mu_phi_ctr, Y$mu_phi_opi, Y$mu_phi_amp,
                        expression(mu[phi]), c(-3, 3))

p_c <- plot_parameter(Y$mu_c_ctr, Y$mu_c_opi, Y$mu_c_amp,
                      expression(mu[c]), c(0, 1))

# Combine plots in a grid
combined_plot <- grid.arrange(
  p_theta, p_delta, p_alpha, p_phi, p_c,
  ncol = 2,
  nrow = 3
)

# Save the plot
ggsave("VSE_group_comparison_parameters_styled.png", 
       combined_plot, 
       width = 10, 
       height = 12, 
       dpi = 300,
       bg = "white")

cat("Saved: VSE_group_comparison_parameters_styled.png\n")

# Print summary statistics
cat("\n=== PARAMETER COMPARISON SUMMARY ===\n")

print_param_summary <- function(param_name, ctr, opi, amp) {
  cat("\n", param_name, ":\n", sep="")
  cat("  Controls:     ", sprintf("%.3f [%.3f, %.3f]", 
                                  median(ctr), 
                                  quantile(ctr, 0.025), 
                                  quantile(ctr, 0.975)), "\n")
  cat("  Heroin:       ", sprintf("%.3f [%.3f, %.3f]", 
                                  median(opi), 
                                  quantile(opi, 0.025), 
                                  quantile(opi, 0.975)), "\n")
  cat("  Amphetamine:  ", sprintf("%.3f [%.3f, %.3f]", 
                                  median(amp), 
                                  quantile(amp, 0.025), 
                                  quantile(amp, 0.975)), "\n")
}

print_param_summary("theta", Y$mu_theta_ctr, Y$mu_theta_opi, Y$mu_theta_amp)
print_param_summary("delta", Y$mu_delta_ctr, Y$mu_delta_opi, Y$mu_delta_amp)
print_param_summary("alpha", Y$mu_alpha_ctr, Y$mu_alpha_opi, Y$mu_alpha_amp)
print_param_summary("phi", Y$mu_phi_ctr, Y$mu_phi_opi, Y$mu_phi_amp)
print_param_summary("c", Y$mu_c_ctr, Y$mu_c_opi, Y$mu_c_amp)


# ============================================================
# Pairwise differences from JAGS (posterior median [95% CI])
# ============================================================

cat("\n=== PAIRWISE DIFFERENCES (POSTERIOR MEDIAN [95% CI]) ===\n")

summarize_posterior <- function(x) {
  
  x <- as.numeric(x)
  
  q <- quantile(x, c(0.025, 0.975))
  
  c(
    median = median(x),
    lCI = q[1],
    uCI = q[2]
  )
}



print_diff <- function(name, x) {
  s <- summarize_posterior(x)
  cat(sprintf("  %-25s %.3f [%.3f, %.3f]\n",
              name, s["median"], s["lCI.2.5%"], s["uCI.97.5%"]))
}

# ---- theta ----
cat("\ntheta:\n")
print_diff("Controls − Heroin:",      Y$diff_theta_ctr_opi)
print_diff("Controls − Amphetamine:", Y$diff_theta_ctr_amp)
print_diff("Heroin − Amphetamine:",   Y$diff_theta_opi_amp)

# ---- delta ----
cat("\ndelta:\n")
print_diff("Controls − Heroin:",      Y$diff_delta_ctr_opi)
print_diff("Controls − Amphetamine:", Y$diff_delta_ctr_amp)
print_diff("Heroin − Amphetamine:",   Y$diff_delta_opi_amp)

# ---- alpha ----
cat("\nalpha:\n")
print_diff("Controls − Heroin:",      Y$diff_alpha_ctr_opi)
print_diff("Controls − Amphetamine:", Y$diff_alpha_ctr_amp)
print_diff("Heroin − Amphetamine:",   Y$diff_alpha_opi_amp)

# ---- phi ----
cat("\nphi:\n")
print_diff("Controls − Heroin:",      Y$diff_phi_ctr_opi)
print_diff("Controls − Amphetamine:", Y$diff_phi_ctr_amp)
print_diff("Heroin − Amphetamine:",   Y$diff_phi_opi_amp)

# ---- c ----
cat("\nc:\n")
print_diff("Controls − Heroin:",      Y$diff_c_ctr_opi)
print_diff("Controls − Amphetamine:", Y$diff_c_ctr_amp)
print_diff("Heroin − Amphetamine:",   Y$diff_c_opi_amp)


# -----------------------------
# Function to plot parameter densities with priors
# -----------------------------
plot_parameter <- function(ctr_samples, opi_samples, amp_samples, 
                           param_label, xlim_range,
                           prior_fun = NULL, prior_args = list()) {
  
  # Posterior data
  df <- data.frame(
    value = c(ctr_samples, opi_samples, amp_samples),
    group = factor(rep(c("Controls", "Heroin", "Amphetamine"), 
                       times = c(length(ctr_samples), length(opi_samples), length(amp_samples))),
                   levels = c("Controls", "Heroin", "Amphetamine"))
  )
  
  # Base posterior density plot
  p <- ggplot(df, aes(x = value, fill = group, color = group)) +
    geom_density(alpha = 0.5, size = 1) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(title = param_label,
         x = "",
         y = "Density") +
    xlim(xlim_range) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
  
  # Add priors if provided
  if(!is.null(prior_fun)) {
    for(group_name in c("Controls", "Heroin", "Amphetamine")) {
      p <- p + stat_function(
        fun = prior_fun,
        args = prior_args[[group_name]],
        color = "grey9",
        linetype = "dashed",
        size = .6
      )
    }
  }
  
  return(p)
}

# -----------------------------
# Define priors (adjust to your JAGS priors)
# -----------------------------
# Beta(1,1) for θ, δ, α, c; Normal(0,1) for φ
prior_beta <- list(Controls = list(shape1=1, shape2=1),
                   Heroin   = list(shape1=1, shape2=1),
                   Amphetamine = list(shape1=1, shape2=1))

prior_normal <- list(Controls = list(mean=0, sd=1),
                     Heroin   = list(mean=0, sd=1),
                     Amphetamine = list(mean=0, sd=1))

# -----------------------------
# Generate individual plots
# -----------------------------
p_theta <- plot_parameter(Y$mu_theta_ctr, Y$mu_theta_opi, Y$mu_theta_amp,
                          expression(mu[theta]), c(0,1),
                          prior_fun = dbeta, prior_args = prior_beta)

p_delta <- plot_parameter(Y$mu_delta_ctr, Y$mu_delta_opi, Y$mu_delta_amp,
                          expression(mu[delta]), c(0,1),
                          prior_fun = dbeta, prior_args = prior_beta)

p_alpha <- plot_parameter(Y$mu_alpha_ctr, Y$mu_alpha_opi, Y$mu_alpha_amp,
                          expression(mu[alpha]), c(0,1),
                          prior_fun = dbeta, prior_args = prior_beta)

p_phi <- plot_parameter(Y$mu_phi_ctr, Y$mu_phi_opi, Y$mu_phi_amp,
                        expression(mu[phi]), c(-3,3),
                        prior_fun = dnorm, prior_args = prior_normal)

p_c <- plot_parameter(Y$mu_c_ctr, Y$mu_c_opi, Y$mu_c_amp,
                      expression(mu[c]), c(0,1),
                      prior_fun = dbeta, prior_args = prior_beta)

# -----------------------------
# Combine plots into a grid
# -----------------------------
combined_plot <- grid.arrange(
  p_theta, p_delta, p_alpha, p_phi, p_c,
  ncol = 2,
  nrow = 3
)

# -----------------------------
# Save figure
# -----------------------------
ggsave("VSE_group_comparison_parameters_with_priors.png",
       combined_plot,
       width = 10,
       height = 12,
       dpi = 300,
       bg = "white")

cat("Saved: VSE_group_comparison_parameters_with_priors.png\n")

