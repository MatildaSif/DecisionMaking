#take argument from bash
arg <- commandArgs(TRUE)
#make it an integer
seed <- as.integer(arg)
#load packs
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm,tidyverse)

set.seed(seed)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 10 # size of our subdivisions for pseudorandomization
freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
bad_r <- 100 # "bad" winnings
bad_freq_l <- -250 # "bad" frequent loss
bad_infreq_l <- -1250 # "bad" infrequent loss
good_r <- 50 # "good" winnings
good_freq_l <- -50 # "good" frequent loss
good_infreq_l <- -250 # "good" infrequent loss

# Bad frequent
A_R <- rep(bad_r, nstruct) # we win on every trials
A_L <- c(rep(bad_freq_l, nstruct*freq),rep(0,nstruct*(1-freq))) # we have losses half of the time

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq))) # we have losses 1/10th of the time

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq),rep(0,nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq)))

# create the pseudorandomized full payoff structure
A <- array(NA,ntrials) # setting up and empty array to be filled
B <- array(NA,ntrials)
C <- array(NA,ntrials)
D <- array(NA,ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}


payoff <- cbind(A,B,C,D)/100 # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with

print(" ---- payoff structure made ---- ")
###--------------Run full parameter recovery -------------
niterations <- 5 # fewer because it takes too long
nsubs <- 48 # mimicking the data structure from Ahn et al.
ntrials_all <- rep(100, 48) # all 48 simulated subs have 100 trials each

# mu
true_mu_a_rew <- array(NA,c(niterations))
true_mu_a_pun <- array(NA,c(niterations))
true_mu_K <- array(NA,c(niterations))
true_mu_theta <- array(NA,c(niterations))
true_mu_omega_f <- array(NA,c(niterations))
true_mu_omega_p <- array(NA,c(niterations))

infer_mu_a_rew <- array(NA,c(niterations))
infer_mu_a_pun <- array(NA,c(niterations))
infer_mu_K <- array(NA,c(niterations))
infer_mu_theta <- array(NA,c(niterations))
infer_mu_omega_f <- array(NA,c(niterations))
infer_mu_omega_p <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_lambda_a_rew <- array(NA,c(niterations))
true_lambda_a_pun <- array(NA,c(niterations))
true_lambda_K <- array(NA,c(niterations))
true_lambda_theta <- array(NA,c(niterations))
true_lambda_omega_f <- array(NA,c(niterations))
true_lambda_omega_p <- array(NA,c(niterations))

infer_lambda_a_rew <- array(NA,c(niterations))
infer_lambda_a_pun <- array(NA,c(niterations))
infer_lambda_K <- array(NA,c(niterations))
infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_omega_f <- array(NA,c(niterations))
infer_lambda_omega_p <- array(NA,c(niterations))


print(" ---- First iteration starts ---- ")
start_time = Sys.time()
print(start_time)
for (i in 1:niterations) {
  ntrials <- ntrials_all
  
  # let's see how robust the model is. Does it recover all sorts of values?
  mu_a_rew <- runif(1,0,1)
  mu_a_pun <- runif(1,0,1)
  mu_K <- runif(1,0,2)
  mu_theta <- runif(1,.2,2) # could also just be a set value (e.g. 1) to simplify the model a bit
  mu_omega_f <- runif(1,-2,2)
  mu_omega_p <- runif(1,-2,2)
  
  sigma_a_rew <- runif(1,0,0.1)
  sigma_a_pun <- runif(1,0,0.1)
  sigma_K <- runif(1,0,0.2)
  sigma_theta <- runif(1,0,0.2) # if theta is just a set value (e.g. 1), then this isn't relevant anymore
  sigma_omega_f <- runif(1,0,0.4)
  sigma_omega_p <- runif(1,0,0.4)
  
  # sigma_a_rew <- runif(1,0,.5)
  # sigma_a_pun <- runif(1,0,.5)
  # sigma_K <- runif(1,0,.5)
  # sigma_theta <- runif(1,0,.5)
  # sigma_omega_f <- runif(1,0,.5)
  # sigma_omega_p <- runif(1,0,.5)
  
  source('hier_ORL_sim.R')
  ORL_sims <- hier_ORL_sim(payoff,nsubs,ntrials,mu_a_rew,mu_a_pun,
                           mu_K,mu_theta,mu_omega_f,mu_omega_p,
                           sigma_a_rew,sigma_a_pun,sigma_K,sigma_theta,
                           sigma_omega_f,sigma_omega_p)
  
  x <- ORL_sims$x
  X <- ORL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","ntrials","nsubs") 
  params<-c("mu_a_rew","mu_a_pun",
            "mu_K","mu_theta","mu_omega_f","mu_omega_p",
            "lambda_a_rew","lambda_a_pun","lambda_K","lambda_theta",
            "lambda_omega_f","lambda_omega_p")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="hier_ORL.txt", n.chains=4, 
                           n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  
  #save output
  if (i ==1) {
    output_df <- tibble(
      #general
      seed = seed,
      n_iteration = i,
      # mu
      true_mu_a_rew = mu_a_rew,
      true_mu_a_pun = mu_a_pun,
      true_mu_K = mu_K,
      true_mu_theta = mu_theta,
      true_mu_omega_f = mu_omega_f,
      true_mu_omega_p = mu_omega_p,
      #inferred mu
      infer_mu_a_rew = MPD(Y$mu_a_rew),
      infer_mu_a_pun = MPD(Y$mu_a_pun),
      infer_mu_K = MPD(Y$mu_K),
      infer_mu_theta = MPD(Y$mu_theta),
      infer_mu_omega_f = MPD(Y$mu_omega_f),
      infer_mu_omega_p = MPD(Y$mu_omega_p),
      # lambda
      true_lambda_a_rew = sigma_a_rew,
      true_lambda_a_pun = sigma_a_pun,
      true_lambda_K = sigma_K,
      true_lambda_theta = sigma_theta,
      true_lambda_omega_f = sigma_omega_f,
      true_lambda_omega_p = sigma_omega_p,
      # find maximum a posteriori
      infer_lambda_a_rew = MPD(Y$lambda_a_rew),
      infer_lambda_a_pun = MPD(Y$lambda_a_pun),
      infer_lambda_K =  MPD(Y$lambda_K),
      infer_lambda_theta = MPD(Y$lambda_theta),
      infer_lambda_omega_f = MPD(Y$lambda_omega_f),
      infer_lambda_omega_p = MPD(Y$lambda_omega_p)
    )
  }else {
    temporary_df <- tibble(
      #general
      seed = seed,
      n_iteration = i,
      # mu
      true_mu_a_rew = mu_a_rew,
      true_mu_a_pun = mu_a_pun,
      true_mu_K = mu_K,
      true_mu_theta = mu_theta,
      true_mu_omega_f = mu_omega_f,
      true_mu_omega_p = mu_omega_p,
      #inferred mu
      infer_mu_a_rew = MPD(Y$mu_a_rew),
      infer_mu_a_pun = MPD(Y$mu_a_pun),
      infer_mu_K = MPD(Y$mu_K),
      infer_mu_theta = MPD(Y$mu_theta),
      infer_mu_omega_f = MPD(Y$mu_omega_f),
      infer_mu_omega_p = MPD(Y$mu_omega_p),
      # lambda
      true_lambda_a_rew = sigma_a_rew,
      true_lambda_a_pun = sigma_a_pun,
      true_lambda_K = sigma_K,
      true_lambda_theta = sigma_theta,
      true_lambda_omega_f = sigma_omega_f,
      true_lambda_omega_p = sigma_omega_p,
      # find maximum a posteriori
      infer_lambda_a_rew = MPD(Y$lambda_a_rew),
      infer_lambda_a_pun = MPD(Y$lambda_a_pun),
      infer_lambda_K =  MPD(Y$lambda_K),
      infer_lambda_theta = MPD(Y$lambda_theta),
      infer_lambda_omega_f = MPD(Y$lambda_omega_f),
      infer_lambda_omega_p = MPD(Y$lambda_omega_p)
    )
    
    output_df <- rbind(output_df,temporary_df)
    
  }
  
  
  
  print(paste("Progress:",i/niterations*100,"%",sep = " "))
  current_time <- Sys.time()
  time_spent <-current_time - start_time
  print(paste("Since start:",time_spent,sep = " "))
  
}


#write output to folder
name <- paste("seed2",seed, sep = "_")
out_f <- paste("./results/",name,".rds", sep = "")
write_rds(output_df,out_f)