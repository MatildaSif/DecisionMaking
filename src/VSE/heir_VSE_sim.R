heir_VSE_sim <- function(R_decks, L_decks, nsubs, ntrials_all, 
                         mu_theta, mu_delta, mu_alpha, mu_phi, mu_c,
                         sigma_theta, sigma_delta, sigma_alpha, sigma_phi, sigma_c) {
  
  
  x <- array(NA, c(nsubs,  ntrials_all[1]))
  R <- array(NA, c(nsubs,  ntrials_all[1]))
  L <- array(NA, c(nsubs,  ntrials_all[1]))
  
  print(dim(R))
  print(dim(L))
  
  for (s in 1:nsubs) {
    # Sample parameters
    theta <- rtruncnorm(1, 0, 1, mu_theta, sigma_theta)      # learning rate
    delta <- rtruncnorm(1, 0, Inf, mu_delta, sigma_delta)    # exploration bonus
    alpha <- rtruncnorm(1, 0, 1, mu_alpha, sigma_alpha)      # exploration weight
    phi <- rtruncnorm(1, 0, 1, mu_phi, sigma_phi)            # decay rate
    c_param <- rtruncnorm(1, 0, Inf, mu_c, sigma_c)                # inverse temp
    
    exploit <- array(0, c(ntrials_all[s], 4))
    explore <- array(0, c( ntrials_all[s], 4))
    Ev <- array(0, c( ntrials_all[s], 4))
    
    # First trial: random choice
    x[s, 1] <- sample(1:4, 1)
    R[s, 1] <- R_decks[1, x[s, 1]]
    L[s, 1] <- L_decks[1, x[s, 1]]
    
    for (t in 2:ntrials_all[s]) {
      
      # Outcome value (can use prospect theory or simple difference) 
      v <- (R[s,t-1]^0.5) - (L[s,t-1]^0.5)
      
      for (d in 1:4) {
        if (d == x[s, t-1]) {
          # CHOSEN: update exploitation, set exploration bonus
          exploit[t, d] <- exploit[t-1, d] + theta * (v - exploit[t-1, d])
          explore[t, d] <- delta
        } else {
          # UNCHOSEN: no exploit update, decay exploration
          exploit[t, d] <- exploit[t-1, d]
          explore[t, d] <- explore[t-1, d] * (1 - phi)
        }
        
        # Combined value
        Ev[t, d] <- exploit[t, d] + alpha * explore[t, d]
      }
      
      # Softmax choice
      exp_Ev <- exp(c_param * Ev[t, ])
      p <- exp_Ev / sum(exp_Ev)
      x[s, t] <- sample(1:4, 1, prob = p)
      
      # Get outcome
      R[s, t] <- R_decks[t, x[s, t]]
      L[s, t] <- L_decks[t, x[s, t]]
    }
  }
  
  return(list(x = x, R = R, L = L))
}