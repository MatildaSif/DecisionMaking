hier_ORL_sim <- function(payoff,nsubs,ntrials,mu_a_rew,mu_a_pun,
                         mu_K,mu_theta,mu_omega_f,mu_omega_p,
                         sigma_a_rew,sigma_a_pun,sigma_K,sigma_theta,
                         sigma_omega_f,sigma_omega_p) {
  
  # arrays to populate for simulation
  x <- array(NA,c(nsubs,ntrials[1]))
  X <- array(NA,c(nsubs,ntrials[1]))
  Ev <- array(NA,c(nsubs,ntrials[1],4))
  
  
  for (s in 1:nsubs) {
    
    # free parameters - sampled from a normal distribution with group mean and sd
    a_rew <- rtruncnorm(1,0,,mu_a_rew,sigma_a_rew)
    a_pun <- rtruncnorm(1,0,,mu_a_pun,sigma_a_pun)
    K <- rtruncnorm(1,0,,mu_K,sigma_K)
    theta <- rtruncnorm(1,0,,mu_theta,sigma_theta)
    omega_f <- rtruncnorm(1,0,,mu_omega_f,sigma_omega_f)
    omega_p <- rtruncnorm(1,0,,mu_omega_p,sigma_omega_p)

    # arrays to populate for simulation 
    Ev_update <- array(NA,c(ntrials[s],4))
    
    signX <- array(NA,c(ntrials[s]))
    Ef_cho <- array(NA,c(ntrials[s],4))
    Ef_not <- array(NA,c(ntrials[s],4))
    Ef <- array(NA,c(ntrials[s],4))
    
    PS <- array(NA,c(ntrials[s],4))
    
    V <- array(NA,c(ntrials[s],4))
    
    exp_p <- array(NA,c(ntrials[s],4))
    p <- array(NA,c(ntrials[s],4))
    
    x[s,1] <- rcat(1,c(.25,.25,.25,.25))
    
    X[s,1] <- payoff[1, x[s,1]]
    
    Ev[s,1,] <- rep(0,4)
    
    Ef[1,] <- rep(0,4)
    
    PS[1,] <- rep(0,4)
    
    for (t in 2:ntrials) {
      
      #this is important mention this as constructing model
      signX[t] <- ifelse(X[s,t-1]<0,-1,1)
      
      for (d in 1:4) {
        
        # -------- Updating expected values ------------------------
        Ev_update[t,d] <- ifelse(X[s,t-1]>=0,
                                 Ev[s,t-1,d] + a_rew*((X[s,t-1]) - Ev[s,t-1,d]), 
                                 Ev[s,t-1,d] + a_pun*((X[s,t-1]) - Ev[s,t-1,d])
        )
        
        Ev[s,t,d] <- ifelse(d==x[s,t-1],Ev_update[t,d],Ev[s,t-1,d])
        
        # -------- Updating expected frequencies ------------------------
        #update expected frequencies for ALL decks - AS IF THEY WERE ALL CHOSEN
        Ef_cho[t,d] <- ifelse(X[s,t-1]>=0, 
                              Ef[t-1,d] + a_rew*(signX[t] - Ef[t-1,d]),
                              Ef[t-1,d] + a_pun*(signX[t] - Ef[t-1,d])
        )
        
        #update expected frequencies for ALL decks - AS IF THEY WERE ALL UNCHOSEN. 
        Ef_not[t,d] <- ifelse(X[s,t-1]>=0, 
                              Ef[t-1,d] + a_pun*(-(signX[t]/3) - Ef[t-1,d]),
                              Ef[t-1,d] + a_rew*(-(signX[t]/3) - Ef[t-1,d])
        ) 
        
        #copy appropriate values to ef variable
        Ef[t,d] <- ifelse(d==x[s,t-1],Ef_cho[t,d],Ef_not[t,d])  
        #------------------------------------------------------------------
        
        #ifelse needed to disctiminate chosen and unchosen decks
        PS[t,d] <- ifelse(x[s,t-1]==d,1/(1+K),PS[t-1,d]/(1+K))
        
        V[t,d] <- Ev[s,t,d] + Ef[t,d]*omega_f + PS[t,d]*omega_p
        
        exp_p[t,d] <- exp(theta*V[t,d])
        
      }
      
      for (d in 1:4) {
        p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
      }
      
    x[s,t] <- rcat(1,p[t,])
      
    X[s,t] <- payoff[t,x[s,t]]
    
    }
  }
  
  result <- list(x=x,
                 X=X,
                 Ev=Ev)
  
  return(result)
  
}


