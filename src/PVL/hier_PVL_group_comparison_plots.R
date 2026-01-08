# Compare PVL posteriors across groups

library(tidyverse)
library(ggplot2)

setwd('/work/JoMat/DecisionMaking/')

output_dir <- "src/PVL/outputs/parameter_estimation"
fig_dir <- file.path(output_dir, "figures")

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}


# Load fitted models
load(file.path(output_dir, "PVL_amphetamine_posteriors.RData"))
Y_amp <- samples$BUGSoutput$sims.list

load(file.path(output_dir, "PVL_heroin_posteriors.RData"))
Y_her <- samples$BUGSoutput$sims.list

load(file.path(output_dir, "PVL_control_posteriors.RData"))
Y_ctl <- samples$BUGSoutput$sims.list


# Combine group-level means
group_means <- bind_rows(
  data.frame(w = Y_amp$mu_w, A = Y_amp$mu_A,
             theta = Y_amp$mu_theta, a = Y_amp$mu_a,
             Group = "Amphetamine"),
  data.frame(w = Y_her$mu_w, A = Y_her$mu_A,
             theta = Y_her$mu_theta, a = Y_her$mu_a,
             Group = "Heroin"),
  data.frame(w = Y_ctl$mu_w, A = Y_ctl$mu_A,
             theta = Y_ctl$mu_theta, a = Y_ctl$mu_a,
             Group = "Control")
)

group_means_long <- pivot_longer(
  group_means,
  cols = c(w, A, theta, a),
  names_to = "Parameter",
  values_to = "Value"
)


# Plot: Group-level μ
# -------------------------------
p_mu <- ggplot(group_means_long,
               aes(x = Value, fill = Group, colour = Group)) +
  geom_density(alpha = 0.35) +
  facet_wrap(~ Parameter, scales = "free", ncol = 2) +
  labs(
    title = "Group-Level Posterior Distributions (μ)",
    x = "Parameter value",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  file.path(fig_dir, "PVL_group_means_overlay.png"),
  p_mu,
  width = 11,
  height = 8,
  dpi = 300
)


# Combine group-level precisions
# -------------------------------
group_prec <- bind_rows(
  data.frame(lambda_w = Y_amp$lambda_w,
             lambda_A = Y_amp$lambda_A,
             lambda_theta = Y_amp$lambda_theta,
             lambda_a = Y_amp$lambda_a,
             Group = "Amphetamine"),
  data.frame(lambda_w = Y_her$lambda_w,
             lambda_A = Y_her$lambda_A,
             lambda_theta = Y_her$lambda_theta,
             lambda_a = Y_her$lambda_a,
             Group = "Heroin"),
  data.frame(lambda_w = Y_ctl$lambda_w,
             lambda_A = Y_ctl$lambda_A,
             lambda_theta = Y_ctl$lambda_theta,
             lambda_a = Y_ctl$lambda_a,
             Group = "Control")
)

group_prec_long <- pivot_longer(
  group_prec,
  cols = starts_with("lambda"),
  names_to = "Parameter",
  values_to = "Value"
)


# Plot: Group-level λ
# -------------------------------
p_lambda <- ggplot(group_prec_long,
                   aes(x = Value, fill = Group, colour = Group)) +
  geom_density(alpha = 0.35) +
  facet_wrap(~ Parameter, scales = "free", ncol = 2) +
  labs(
    title = "Group-Level Posterior Distributions (λ)",
    x = "Precision",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  file.path(fig_dir, "PVL_group_precisions_overlay.png"),
  p_lambda,
  width = 11,
  height = 8,
  dpi = 300
)


# Subject-level posterior means
# -------------------------------
subject_means <- bind_rows(
  data.frame(w = colMeans(Y_amp$w),
             A = colMeans(Y_amp$A),
             theta = colMeans(Y_amp$theta),
             a = colMeans(Y_amp$a),
             Group = "Amphetamine"),
  data.frame(w = colMeans(Y_her$w),
             A = colMeans(Y_her$A),
             theta = colMeans(Y_her$theta),
             a = colMeans(Y_her$a),
             Group = "Heroin"),
  data.frame(w = colMeans(Y_ctl$w),
             A = colMeans(Y_ctl$A),
             theta = colMeans(Y_ctl$theta),
             a = colMeans(Y_ctl$a),
             Group = "Control")
)

subject_means_long <- pivot_longer(
  subject_means,
  cols = c(w, A, theta, a),
  names_to = "Parameter",
  values_to = "Value"
)

p_subjects <- ggplot(subject_means_long,
                     aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Parameter, scales = "free") +
  labs(
    title = "Distribution of Subject-Level Posterior Means",
    x = "Posterior mean",
    y = "Density"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(fig_dir, "PVL_subject_level_means_overlay.png"),
  p_subjects,
  width = 11,
  height = 8,
  dpi = 300
)

cat("All comparison plots saved to:", fig_dir, "\n")
