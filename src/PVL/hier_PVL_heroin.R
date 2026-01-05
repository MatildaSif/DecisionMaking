install.packages("pacman")
pacman::p_load(R2jags, parallel, ggplot2)

set.seed(1983)
setwd('/work/JoMat/DecisionMaking/')

# Create output folder structure in PVL directory
output_dir <- "PVL/outputs/parameter_estimation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# defining a function for calculating the maximum of the posterior density
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
ctr_data <- read.table("data/IGTdata_heroin.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# identify and count unique subject IDs
subIDs <- unique(ctr_data$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# all choices (x) and outcomes (X)
x_raw <- ctr_data$deck
X_raw <- ctr_data$gain + ctr_data$loss #note the sign!

#--- assign choices and outcomes in trial x sub matrix
# empty arrays to fill
ntrials_all <- array(0,c(nsubs))
x_all <- array(0,c(nsubs,ntrials_max))
X_all <- array(0,c(nsubs,ntrials_max))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  ntrials_all[s] <- length(x_raw[ctr_data$subjID==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[ctr_data$subjID==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[ctr_data$subjID==subIDs[s]] 
  length(X_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  X_all[s,] <- X_sub
  
}

###########################################################
#---------- run the hierarchical model on controls --------
###########################################################
x <- x_all
X <- X_all
ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 

# MODIFIED: Include subject-level parameters AND log_lik for model comparison
params <- c("mu_w", "mu_A", "mu_theta", "mu_a",
            "lambda_w", "lambda_A", "lambda_theta", "lambda_a",
            "w", "A", "theta", "a")                 

cat("Fitting hierarchical PVL model to heroin addicts... \n")
cat("Number of subjects:", nsubs, "\n")

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="src/PVL/hier_PVL.txt",
                         n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
end_time = Sys.time()
cat("Fitting time:", end_time - start_time, "\n")

# Check dimensions
cat("\n=== Checking saved parameters ===\n")
Y <- samples$BUGSoutput$sims.list
cat("Group-level parameters:\n")
cat("  mu_w:", length(Y$mu_w), "samples\n")
cat("  mu_A:", length(Y$mu_A), "samples\n")
cat("  mu_theta:", length(Y$mu_theta), "samples\n")
cat("  mu_a:", length(Y$mu_a), "samples\n")

cat("\nSubject-level parameters:\n")
cat("  w:", dim(Y$w), "(samples x subjects)\n")
cat("  A:", dim(Y$A), "(samples x subjects)\n")
cat("  theta:", dim(Y$theta), "(samples x subjects)\n")
cat("  a:", dim(Y$a), "(samples x subjects)\n")

if ("log_lik" %in% names(Y)) {
  cat("\nLog-likelihood:", dim(Y$log_lik), "(samples x subjects)\n")
}

# Save group-level posterior plots
png(file.path(output_dir, "01_group_level_posteriors_heroin.png"), width=800, height=800)
par(mfrow=c(2,2))
plot(density(Y$mu_w), main="Group mean: w")
plot(density(Y$mu_A), main="Group mean: A")
plot(density(Y$mu_theta), main="Group mean: theta")
plot(density(Y$mu_a), main="Group mean: a")
dev.off()
cat("Saved: 01_group_level_posteriors_heroin.png\n")

# Save group-level precision plots
png(file.path(output_dir, "02_group_level_precisions_heroin.png"), width=800, height=800)
par(mfrow=c(2,2))
plot(density(Y$lambda_w), main="Group precision: lambda_w")
plot(density(Y$lambda_A), main="Group precision: lambda_A")
plot(density(Y$lambda_theta), main="Group precision: lambda_theta")
plot(density(Y$lambda_a), main="Group precision: lambda_a")
dev.off()
cat("Saved: 02_group_level_precisions_heroin.png\n")

# Save subject-level posterior plots (Subject 1)
png(file.path(output_dir, "03_subject_1_posteriors_heroin.png"), width=800, height=800)
par(mfrow=c(2,2))
plot(density(Y$w[,1]), main="Subject 1: w")
plot(density(Y$A[,1]), main="Subject 1: A")
plot(density(Y$theta[,1]), main="Subject 1: theta")
plot(density(Y$a[,1]), main="Subject 1: a")
dev.off()
cat("Saved: 03_subject_1_posteriors_heroin.png\n")

# Save everything to RData file
rdata_file <- file.path(output_dir, "PVL_heroin_posteriors.RData")
save(samples, x_all, X_all, ntrials_all, nsubs, subIDs,
     file = rdata_file) 
cat("\n=== Saved posteriors to", rdata_file, "===\n")

# Save summary statistics to text file
summary_file <- file.path(output_dir, "model_summary_heroin.txt")
sink(summary_file)

cat("=== FITTED PARAMETERS (Heroin) ===\n")
cat("\nGroup Means (mu):\n")
cat("w:    ", round(MPD(Y$mu_w), 3), "\n")
cat("A:    ", round(MPD(Y$mu_A), 3), "\n")
cat("theta:", round(MPD(Y$mu_theta), 3), "\n")
cat("a:    ", round(MPD(Y$mu_a), 3), "\n")

cat("\nGroup Precisions (lambda):\n")
cat("w:    ", round(MPD(Y$lambda_w), 3), "\n")
cat("A:    ", round(MPD(Y$lambda_A), 3), "\n")
cat("theta:", round(MPD(Y$lambda_theta), 3), "\n")
cat("a:    ", round(MPD(Y$lambda_a), 3), "\n")

cat("\nGroup SDs (converted from lambda):\n")
cat("w:    ", round(MPD(1/sqrt(Y$lambda_w)), 3), "\n")
cat("A:    ", round(MPD(1/sqrt(Y$lambda_A)), 3), "\n")
cat("theta:", round(MPD(1/sqrt(Y$lambda_theta)), 3), "\n")
cat("a:    ", round(MPD(1/sqrt(Y$lambda_a)), 3), "\n")

cat("\n=== MODEL FIT ===\n")
cat("DIC:", round(samples$BUGSoutput$DIC, 2), "\n")
cat("pD: ", round(samples$BUGSoutput$pD, 2), "\n")

sink()
cat("Saved: model_summary_heroin.txt\n")
