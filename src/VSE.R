VSE <- function(R,L,ntrials,theta,delta,alpha,phi,c) {

  # arrays to populate for simulation
  x <- array(0,c(ntrials))
  X <- array(0,c(ntrials))

  v <- array(0,c(ntrials))
  exploit <- array(0,c(ntrials,4))
  explore <- array(0,c(ntrials,4))

  exp_p <- array(0,c(ntrials,4))
  p <- array(0,c(ntrials,4))
  
  # free parameters - turn back on when constructing
  # theta <- .5
  # delta <- .5
  # alpha <- .3
  # phi <- .5
  # c <- 1
  
  # need to initialise first choice in this model. Add this!!!
  x[1] <- sample(c(1:4),1)
  
  for (t in 2:ntrials) {
    
    # equation 1 - no prospect theory - discuss
    # efficient code, R and L are indexed by choice on last trial
    v[t] <- (R[t,x[t-1]]^theta) - (L[t,x[t-1]]^theta)
    
    for (d in 1:4) {
      
      #---equation 2--- update exploitation. 
      exploit[t,d] <- ifelse(d==x[t-1],
                             exploit[t-1,d]*delta+v[t],
                             exploit[t-1,d]*delta)
      #---equation 3--- 
      explore[t,d] <- ifelse(d==x[t-1],
                             0,
                             explore[t-1,d] + (alpha*(phi-explore[t-1,d])))
      
      #---equation 4--- 
      exp_p[t,d] <- exp((explore[t,d]+exploit[t,d])*c)
      
    }
    
    for (d in 1:4) {
      p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
    }
      
    x[t] <- rcat(1,p[t,])
    
  }
  
  result <- list(x=x,
                 X=X,
                 explore=explore,
                 exploit=exploit)
  
  return(result)
  
  
  #turn back on when building
  par(mfrow=c(2,2))
  plot(explore[,1])
  plot(explore[,2])
  plot(explore[,3])
  plot(explore[,4])
  
  plot(exploit[,1])
  plot(exploit[,2])
  plot(exploit[,3])
  plot(exploit[,4])
  plot(x)
  
}
