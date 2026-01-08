library(R2jags)
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("/work/JoMat/DecisionMaking/src/VSE")

# ============================
# Output folder
# ============================
output_dir <- "Convergence"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================
# Load RData
# ============================
load("VSE_group_comparison_results.RData")  # contains 'samples'

# ============================
# Define groups
# ============================
PARAM_MATCH <- c("ctr", "opi", "amp")                # parameter suffixes
FILE_NAME <- c("ctr", "opi", "amphetamine")         # file-friendly names
GROUP_TITLES <- c("Healthy Controls", "Heroin Users", "Amphetamine Users")

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

# ---- Print DIC (fixed) ----
cat("Deviance Information Criterion (DIC):", round(samples$BUGSoutput$DIC, 2), "\n\n")

# ---- Prepare summary table ----
group_params <- grep("^mu_", param_names, value = TRUE)

summary_table <- data.frame(
  Parameter = character(),
  Group = character(),
  Mean = numeric(),
  SD = numeric(),
  Rhat = numeric(),
  stringsAsFactors = FALSE
)

for (g in seq_along(PARAM_MATCH)) {
  param_suffix <- PARAM_MATCH[g]
  group_title <- GROUP_TITLES[g]          # display name
  file_name <- FILE_NAME[g]               # for file outputs
  
  # Select parameters for this group
  group_params_g <- grep(paste0("_", param_suffix, "$"), group_params, value = TRUE)
  
  if(length(group_params_g) > 0){
    s <- samples$BUGSoutput$summary[group_params_g, ]
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
for (g in seq_along(PARAM_MATCH)) {
  param_suffix <- PARAM_MATCH[g]
  group_title <- GROUP_TITLES[g]
  file_name <- FILE_NAME[g]
  
  # Select parameters for this group
  group_params_g <- grep(paste0("_", param_suffix, "$"), group_params, value = TRUE)
  
  if(length(group_params_g) > 0){
    trace_plots <- list()
    for(param in group_params_g){
      trace_plots[[param]] <- create_trace_plot(mcmc_array[, , param], param, rhat[param])
    }
    
    p_traces <- plot_grid(plotlist = trace_plots, ncol = 2)
    trace_file <- file.path(output_dir, paste0("trace_plots_", file_name, ".png"))
    ggsave(trace_file, plot = p_traces, width = 14, height = 10, dpi = 300)
    cat("\n✓ Trace plots saved for", group_title, ":", trace_file, "\n")
  }
}
