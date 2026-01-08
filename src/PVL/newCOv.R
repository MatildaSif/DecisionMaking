library(R2jags)
library(ggplot2)
library(cowplot)
library(dplyr)

# ============================
# Set working directory & output
# ============================
setwd("/work/JoMat/DecisionMaking/")
output_dir <- "src/PVL/outputs/group_comparison/convergence_diagnostics"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================
# Load PVL group comparison results
# ============================
load("src/PVL/outputs/group_comparison/PVL_group_comparison_results.RData")  # contains 'samples'

# ============================
# Define groups
# ============================
PARAM_MATCH <- list(
  control = c("mu_w_ctr", "mu_A_ctr", "mu_theta_ctr", "mu_a_ctr",
              "lambda_w_ctr", "lambda_A_ctr", "lambda_theta_ctr", "lambda_a_ctr"),
  heroin = c("mu_w_her", "mu_A_her", "mu_theta_her", "mu_a_her",
             "lambda_w_her", "lambda_A_her", "lambda_theta_her", "lambda_a_her"),
  amphetamine = c("mu_w_amp", "mu_A_amp", "mu_theta_amp", "mu_a_amp",
                  "lambda_w_amp", "lambda_A_amp", "lambda_theta_amp", "lambda_a_amp")
)

GROUP_TITLES <- c("Healthy Controls", "Heroin Users", "Amphetamine Users")
GROUP_KEYS <- names(PARAM_MATCH)

# ============================
# Trace plot function
# ============================
create_trace_plot <- function(param_array, param_name, rhat_val) {
  n_iter <- nrow(param_array)
  n_chains <- ncol(param_array)
  
  trace_df <- data.frame(
    iteration = rep(1:n_iter, n_chains),
    value = as.vector(param_array),
    chain = factor(rep(1:n_chains, each = n_iter))
  )
  
  chain_colors <- scales::hue_pal()(n_chains)
  
  ggplot(trace_df, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    labs(title = paste("Trace Plot:", param_name),
         subtitle = paste("R-hat =", round(rhat_val, 4)),
         x = "Iteration", y = "Parameter Value") +
    scale_color_manual(values = setNames(chain_colors, 1:n_chains)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ============================
# MODEL SUMMARY
# ============================
rhat <- samples$BUGSoutput$summary[, "Rhat"]
mcmc_array <- samples$BUGSoutput$sims.array
param_names <- rownames(samples$BUGSoutput$summary)

cat("Deviance Information Criterion (DIC):", round(samples$BUGSoutput$DIC, 2), "\n\n")

# ============================
# Prepare summary table per group
# ============================
summary_table <- data.frame(
  Parameter = character(),
  Group = character(),
  Mean = numeric(),
  SD = numeric(),
  Rhat = numeric(),
  stringsAsFactors = FALSE
)

for (g in seq_along(GROUP_KEYS)) {
  group_key <- GROUP_KEYS[g]
  group_title <- GROUP_TITLES[g]
  params_g <- PARAM_MATCH[[group_key]]
  
  # Select parameters that exist in summary
  params_g <- params_g[params_g %in% param_names]
  
  if(length(params_g) > 0){
    s <- samples$BUGSoutput$summary[params_g, ]
    summary_table <- rbind(
      summary_table,
      data.frame(
        Parameter = rownames(s),
        Group = group_title,
        Mean = round(s[, "mean"], 3),
        SD = round(s[, "sd"], 3),
        Rhat = round(s[, "Rhat"], 3),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("=== GROUP-LEVEL PARAMETER SUMMARY ===\n")
print(summary_table)

high_rhat <- summary_table %>% filter(Rhat > 1.1)
if(nrow(high_rhat) > 0){
  cat("\n⚠ Parameters with R-hat > 1.1:\n")
  print(high_rhat)
} else {
  cat("\n✓ All R-hat values <= 1.1\n")
}

# ============================
# TRACE PLOTS (one file per group)
# ============================
for (g in seq_along(GROUP_KEYS)) {
  group_key <- GROUP_KEYS[g]
  group_title <- GROUP_TITLES[g]
  params_g <- PARAM_MATCH[[group_key]]
  params_g <- params_g[params_g %in% dimnames(mcmc_array)[[3]]]
  
  if(length(params_g) > 0){
    trace_plots <- list()
    for(param in params_g){
      trace_plots[[param]] <- create_trace_plot(mcmc_array[, , param], param, rhat[param])
    }
    
    p_traces <- plot_grid(plotlist = trace_plots, ncol = 2)
    trace_file <- file.path(output_dir, paste0("trace_plots_", group_key, ".png"))
    ggsave(trace_file, plot = p_traces, width = 14, height = 10, dpi = 300)
    cat("\n✓ Trace plots saved for", group_title, ":", trace_file, "\n")
  }
}

cat("\n✅ PVL Convergence Diagnostics Completed\nOutputs saved in:", output_dir, "\n")
