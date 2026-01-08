heir_VSE_sim <- function(payoff, nsubs, ntrials,
                         mu_theta, mu_delta, mu_alpha, mu_phi, mu_c,
                         sigma_theta, sigma_delta, sigma_alpha, sigma_phi, sigma_c) {
  
  # Arrays to fill
  x <- array(NA, c(nsubs, max(ntrials)))
  X <- array(NA, c(nsubs, max(ntrials)))
  
  # Store individual parameters
  theta_individual <- array(NA, nsubs)
  delta_individual <- array(NA, nsubs)
  alpha_individual <- array(NA, nsubs)
  phi_individual <- array(NA, nsubs)
  c_individual <- array(NA, nsubs)
  
  for (s in 1:nsubs) {
    # Sample parameters
    theta <- rtruncnorm(1, 0, 1, mu_theta, sigma_theta)
    delta <- rtruncnorm(1, 0, 1, mu_delta, sigma_delta)
    alpha <- rtruncnorm(1, 0, 1, mu_alpha, sigma_alpha)
    phi <- rnorm(1, mu_phi, sigma_phi)
    c <- rtruncnorm(1, 0, 5, mu_c, sigma_c)
    
    # Store individual parameters
    theta_individual[s] <- theta
    delta_individual[s] <- delta
    alpha_individual[s] <- alpha
    phi_individual[s] <- phi
    c_individual[s] <- c
    
    # Initialize arrays
    exploit <- array(0, c(ntrials[s], 4))
    explore <- array(0, c(ntrials[s], 4))
    Ev <- array(0, c(ntrials[s], 4))
    
    # First trial: random choice
    x[s, 1] <- sample(1:4, 1)
    X[s, 1] <- payoff[1, x[s, 1]]
    
    for (t in 2:ntrials[s]) {
      # Get outcome from previous trial
      outcome <- X[s, t-1]
      
      # Separate into reward and loss components
      R <- max(0, outcome)
      L <- abs(min(0, outcome))
      
      # v = R^θ - L^θ
      v <- R^theta - L^theta
      
      for (d in 1:4) {
        if (d == x[s, t-1]) {
          # CHOSEN DECK
          exploit[t, d] <- delta * exploit[t-1, d] + v
          explore[t, d] <- 0
        } else {
          # UNCHOSEN DECK
          exploit[t, d] <- delta * exploit[t-1, d]
          explore[t, d] <- explore[t-1, d] + alpha * (phi - explore[t-1, d])
        }
        
        # Combined value
        Ev[t, d] <- exploit[t, d] + explore[t, d]
      }
      
      # Numerically stable softmax
      scaled_Ev <- c * Ev[t, ]
      max_Ev <- max(scaled_Ev)
      exp_Ev <- exp(scaled_Ev - max_Ev)
      
      if (any(is.na(exp_Ev)) || any(is.infinite(exp_Ev)) || sum(exp_Ev) == 0) {
        p <- rep(0.25, 4)
      } else {
        p <- exp_Ev / sum(exp_Ev)
      }
      
      # Make choice and get payoff
      x[s, t] <- sample(1:4, 1, prob = p)
      X[s, t] <- payoff[t, x[s, t]]
    }
  }
  
  return(list(
    x = x, 
    X = X,
    theta_individual = theta_individual,
    delta_individual = delta_individual,
    alpha_individual = alpha_individual,
    phi_individual = phi_individual,
    c_individual = c_individual
  ))
}