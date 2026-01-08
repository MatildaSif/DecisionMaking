# Fit hierarchical PVL model for all groups 


library(pacman)
pacman::p_load(R2jags, parallel)

set.seed(1983)
setwd('/work/JoMat/DecisionMaking/')


# Group definitions
groups <- list(
  amphetamine = "data/IGTdata_amphetamine_clean.txt",
  heroin      = "data/IGTdata_heroin_clean.txt",
  control     = "data/IGTdata_healthy_control_clean.txt"
)


# Output directory
output_dir <- "src/PVL/outputs/parameter_estimation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Helper function: prepare data
prepare_jags_data <- function(data_file) {
  
  d <- read.table(data_file, header = TRUE)
  
  subIDs <- unique(d$subjID)
  nsubs  <- length(subIDs)
  ntrials_max <- 100
  
  x_raw <- d$deck
  X_raw <- d$gain + d$loss
  
  ntrials_all <- numeric(nsubs)
  x_all <- array(NA, c(nsubs, ntrials_max))
  X_all <- array(NA, c(nsubs, ntrials_max))
  
  for (s in seq_len(nsubs)) {
    
    idx <- d$subjID == subIDs[s]
    ntrials_all[s] <- sum(idx)
    
    x_sub <- x_raw[idx]
    X_sub <- X_raw[idx]
    
    length(x_sub) <- ntrials_max
    length(X_sub) <- ntrials_max
    
    x_all[s, ] <- x_sub
    X_all[s, ] <- X_sub
  }
  
  list(
    x = x_all,
    X = X_all,
    ntrials = ntrials_all,
    nsubs = nsubs,
    subIDs = subIDs
  )
}


# Parameters to monitor
params <- c(
  "mu_w", "mu_A", "mu_theta", "mu_a",
  "lambda_w", "lambda_A", "lambda_theta", "lambda_a",
  "w", "A", "theta", "a", "p"
)


# Loop over groups
for (grp in names(groups)) {
  
  cat("\n====================================\n")
  cat("Fitting PVL model for group:", grp, "\n")
  cat("====================================\n")
  
  # Prepare data
  data_prepared <- prepare_jags_data(groups[[grp]])
  
  data_jags <- list(
    x = data_prepared$x,
    X = data_prepared$X,
    ntrials = data_prepared$ntrials,
    nsubs = data_prepared$nsubs
  )
  
  start_time <- Sys.time()
  
  samples <- jags.parallel(
    data_jags,
    inits = NULL,
    params,
    model.file = "src/PVL/original/hier_PVL.txt",
    n.chains = 4,
    n.iter = 5000,
    n.burnin = 2000,
    n.thin = 1,
    n.cluster = 4
  )
  
  end_time <- Sys.time()
  
  cat("Fitting time:", end_time - start_time, "\n")
  
  # Save output
  save(
    samples,
    data_prepared,
    file = file.path(
      output_dir,
      paste0("PVL_", grp, "_posteriors.RData")
    )
  )
  
  cat("Saved: PVL_", grp, "_posteriors.RData\n", sep = "")
}

cat("\n✅ All group models fitted successfully.\n")
