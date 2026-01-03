
install.packages("pacman")
pacman::p_load(extraDistr, R2jags)


setwd("/work/JoMat/DecisionMaking/src")
set.seed(1702)

#------ create task environment -------------------
#Bad frequent
A_R <- rep(100,10)
A_L <- c(rep(-250,5),rep(0,5))
#Bad infrequent
B_R <- rep(100,10)
B_L <- c(rep(-1250,1),rep(0,9))
#Good frequent
C_R <- rep(50,10)
C_L <- c(rep(-50,5),rep(0,5))
#Good infrequent
D_R <- rep(50,10)
D_L <- c(rep(-250,1),rep(0,9))

A <- c()
for (i in 1:10) {A <- (append(A,A_R + sample(A_L)))}
B <- c()
for (i in 1:10) {B <- (append(B,B_R + sample(B_L)))}
C <- c()
for (i in 1:10) {C <- (append(C,C_R + sample(C_L)))}
D <- c()
for (i in 1:10) {D <- (append(D,D_R + sample(D_L)))}

# NOTE THAT THIS IS DIFFERENT TO OTHER MODELS. WE DON'T HAVE A SINGLE
# PAYOFF SCHEDULE!!!! WE HAVE SEPERATE R AND L SCHEDULES!!!!

# extract rewards from reward trials
RA <- A
RA[A < 0] = 0
RB <- B
RB[B < 0] = 0
RC <- C
RC[C < 0] = 0
RD <- D
RD[D < 0] = 0

# extract losses from reward trials
LA <- A
LA[A > 0] = 0
LB <- B
LB[B > 0] = 0
LC <- C
LC[C > 0] = 0
LD <- D
LD[D > 0] = 0

R <- cbind(RA,RB,RC,RD)
L <- cbind(LA,LB,LC,LD)
L <- abs(L)
R <- R
# THIS IS NOW WHAT WE'LL NEED TO INPUT INTO OUR VSE FUNCTION

ntrials <- 100

#----------------------------------------------------


#-------test VSE function and jags script ---------

#---set params

theta <- .5
delta <- .5
alpha <- .3
phi <- .5
c <- 1

source("VSE.R")
VSE_sims <- VSE(R,L,ntrials,theta,delta,alpha,phi,c)

par(mfrow=c(2,2))
plot(VSE_sims$explore[,1])
plot(VSE_sims$explore[,2])
plot(VSE_sims$explore[,3])
plot(VSE_sims$explore[,4])

plot(VSE_sims$exploit[,1])
plot(VSE_sims$exploit[,2])
plot(VSE_sims$exploit[,3])
plot(VSE_sims$exploit[,4])
plot(VSE_sims$x)

x <- VSE_sims$x
# we don't need payoffs anymore!

# set up jags and run jags model
data <- list("R","L","x","ntrials")  #need R and L in jags model though!
params<-c("theta","delta","alpha","phi","c")
samples <- jags(data, inits=NULL, params,
                model.file ="VSE.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)


###--------------Run full parameter recovery -------------
niterations <- 100 # fewer because it takes too long

true_theta <- array(0,c(niterations))
true_delta <- array(0,c(niterations))
true_alpha <- array(0,c(niterations))
true_phi <- array(0,c(niterations))
true_c <- array(0,c(niterations))

infer_theta <- array(0,c(niterations))
infer_delta <- array(0,c(niterations))
infer_alpha <- array(0,c(niterations))
infer_phi <- array(0,c(niterations))
infer_c <- array(0,c(niterations))

for (i in 1:niterations) {
  
  # # let's see how robust the model is. Does it recover all sorts of values?
  # theta <- runif(1,0,.8) #NOTE IF BOTH LEARNING RATES ARE TOO HIGH, CHOICE RULE WILL CRASH!!!
  # delta <- runif(1,0,.8)
  # alpha <- runif(1,0,1)
  # phi <- rnorm(1,0,.01)
  # c <- runif(1,.5,3)
  
  # Sample from realistic ranges based on what you've observed
  theta <- runif(1, 0.3, 0.7)      # decay: observed range ~0.4-0.6
  delta <- runif(1, 0.4, 1.2)      # exploration bonus: ~0.5-1.0 
  alpha <- runif(1, 0.2, 0.5)      # sensitivity: ~0.25-0.35
  c <- runif(1, 0.6, 1.1)          # consistency: ~0.7-0.95 (NOT 0.5-3!)
  phi <- runif(1, 0.3, 0.6)        # exploration update: ~0.4-0.5
  
  # Truncate to keep in valid range
  theta <- pmax(0.01, pmin(0.99, theta))
  delta <- pmax(0.01, pmin(0.99, delta))

  VSE_sims <- VSE(R,L,ntrials,theta,delta,alpha,phi,c)
  
  x <- VSE_sims$x

  # set up jags and run jags model
  data <- list("R","L","x","ntrials")  #need R and L in jags model though!
  params<-c("theta","delta","alpha","phi","c")
  samples <- jags(data, inits=NULL, params,
                  model.file ="VSE.txt",
                  n.chains=3, n.iter=3000, n.burnin=1000, n.thin=1)
  
  
  true_theta[i] <- theta
  true_delta[i] <- delta
  true_alpha[i] <- alpha
  true_phi[i] <- phi
  true_c[i] <- c
  
  # find maximum a posteriori
  X <- samples$BUGSoutput$sims.list$theta
  infer_theta[i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]

  X <- samples$BUGSoutput$sims.list$delta
  infer_delta[i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
  X <- samples$BUGSoutput$sims.list$alpha
  infer_alpha[i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
  X <- samples$BUGSoutput$sims.list$phi
  infer_phi[i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
  X <- samples$BUGSoutput$sims.list$c
  infer_c[i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
  print(i)
  
}

# let's look at some scatter plots

# Save all plots in one PDF file
pdf("parameter_recovery_plots_VSE_single.pdf", width=10, height=8)
par(mfrow=c(2,2))
plot(true_theta, infer_theta, main="Theta Recovery")
plot(true_delta, infer_delta, main="Delta Recovery")
plot(true_alpha, infer_alpha, main="Alpha Recovery")
plot(true_phi, infer_phi, main="Phi Recovery")
plot(true_c, infer_c, main="C Recovery")
dev.off()




