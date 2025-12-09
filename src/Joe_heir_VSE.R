#Joe's version trying to add all data at once

library(R2jags)

setwd("/work/JohanneSejrskildRejsenhus#9686/DecisionMaking/data")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
ctr_data <- read.table("IGTdata_healthy_control.txt",header=TRUE)
opi_data <- read.table("IGTdata_heorin.txt",header=TRUE)
amp_data <- read.table("IGTdata_amphetamine.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice, gain, and loss ----
# identify and count unique subject IDs
subIDs <- unique(ctr_data$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# all choices (x) and outcomes (X)
x_raw <- ctr_data$deck
R_raw <- ctr_data$gain  
L_raw <- abs(ctr_data$loss) 

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_all <- array(0,c(nsubs))
x_all <- array(0,c(nsubs,ntrials_max))
R_all <- array(0,c(nsubs,ntrials_max))
L_all <- array(0,c(nsubs,ntrials_max))

for (s in 1:nsubs) {

  #record n trials for subject s
  ntrials_all[s] <- length(x_raw[ctr_data$subjID==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[ctr_data$subjID==subIDs[s]] 
  length(x_sub) <- ntrials_max

  R_sub <- R_raw[ctr_data$subjID==subIDs[s]] 
  length(R_sub) <- ntrials_max
  
  L_sub <- L_raw[ctr_data$subjID==subIDs[s]] 
  length(L_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  R_all[s,] <- R_sub
  L_all[s,] <- L_sub
  
}

# Now we'll fit one subject just to make sure everything works

setwd("C:/Users/au199986/Dropbox/Courses/F20/CognitiveModeling/Module4")

x <- x_all[1,]
R <- R_all[1,]
L <- L_all[1,]

ntrials <- ntrials_all[1]

# set up jags and run jags model on one subject
data <- list("x","R","L","ntrials") 
params<-c("theta","delta","alpha","phi","c","p")
samples <- jags(data, inits=NULL, params,
                model.file ="VSE_data.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

# let's look at the posteriors for the parameters
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$theta))
plot(density(samples$BUGSoutput$sims.list$delta))
plot(density(samples$BUGSoutput$sims.list$alpha))
plot(density(samples$BUGSoutput$sims.list$phi))
plot(density(samples$BUGSoutput$sims.list$c))


# Question: how would you expect the data to look on the basis of these posteriors?

# Posterior prediction - start by looking at posteriors for p parameter

p_post <- samples$BUGSoutput$sims.list$p

#plot probability of each deck on trial 32
par(mfrow=c(2,2))
plot(density(p_post[,32,1]))
plot(density(p_post[,32,2]))
plot(density(p_post[,32,3]))
plot(density(p_post[,32,4]))

# which option will be chosen?
x[32]

# is this a good prediction?

# let's write a loop that loop and see how the model goes at predicting responses for all trials 
x_predict <- array(c(ntrials))
for (t in 1:ntrials) {
  
  p_predict <- c(
    density(p_post[,t,1])$x[which(density(p_post[,t,1])$y==max(density(p_post[,t,1])$y))],
    density(p_post[,t,2])$x[which(density(p_post[,t,2])$y==max(density(p_post[,t,2])$y))],
    density(p_post[,t,3])$x[which(density(p_post[,t,3])$y==max(density(p_post[,t,3])$y))],
    density(p_post[,t,4])$x[which(density(p_post[,t,4])$y==max(density(p_post[,t,4])$y))])
  
  x_predict[t] <- which.max(p_predict)
  
}
sum(x_predict==x)

# let's see how the model goes for more than 1 subject. You should run this on all subjects
pred_success <- array(c(nsubs))
for (s in 1:nsubs) {
  
  # fit jags model. Cut and paste from above. Change 1s to s for data and run all.
  x <- x_all[s,]
  R <- R_all[s,]
  L <- L_all[s,]
  
  ntrials <- ntrials_all[s]
  
  # set up jags and run jags model on one subject
  data <- list("x","R","L","ntrials") 
  params<-c("theta","delta","alpha","phi","c","p")
  samples <- jags(data, inits=NULL, params,
                  model.file ="VSE_data.txt",
                  n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
  p_post <- samples$BUGSoutput$sims.list$p
  
  
  
  x_predict <- array(c(ntrials))
  for (t in 1:ntrials) {
  
    p_predict <- c(
        density(p_post[,t,1])$x[which(density(p_post[,t,1])$y==max(density(p_post[,t,1])$y))],
        density(p_post[,t,2])$x[which(density(p_post[,t,2])$y==max(density(p_post[,t,2])$y))],
        density(p_post[,t,3])$x[which(density(p_post[,t,3])$y==max(density(p_post[,t,3])$y))],
        density(p_post[,t,4])$x[which(density(p_post[,t,4])$y==max(density(p_post[,t,4])$y))])
  
    x_predict[t] <- which.max(p_predict)
  
  }
  
  # how many trials did the model predict correctly?
  pred_success[s] <- sum(x_predict==x,na.rm=TRUE)
  print(s)
}



###########################################################
#---------- run the hierarchical model on controls --------
###########################################################
setwd("C:/Users/au199986/Dropbox/Courses/F20/CognitiveModeling/Module4")

x <- x_all
L <- L_all
R <- R_all

ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","L","R","ntrials","nsubs") 
params<-c("mu_theta","mu_delta","mu_alpha","mu_phi","mu_c",
          "lambda_theta","lambda_delta","lambda_alpha","lambda_phi","lambda_c")

samples <- jags(data, inits=NULL, params,
                model.file ="heir_VSE.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

# we'll need to crash out of this. Show slides with posteriors and DIC.

