#### Model simulations: ####
# simulate response and outcomes anew

# ============================================================================================================================= #
#### modSim_M01: DDM (chance performance with bias) ####

modSim_M01 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA # will stay NA
  data$QNoGo <- NA # will stay NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    for (iTrial in 1:nTrial){ # iTrial <- 1
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      
      ## Simulate response and RT using rwiener:
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = par[iSub, 3], delta = par[iSub, 4])
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow responses for updating
      
      # ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      # if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
    } # end iTrial
  } # end iSub
  
  return(data)

}

# ============================================================================================================================= #
#### modSim_02: RL-DDM (learning): #### 

modSim_M02 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){ # iTrial <- 1
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta = par[iSub, 4] + (Q[s, 1] - Q[s, 2]) * par[iSub, 5]
      
      ## Simulate response and RT using rwiener:
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save; crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 6] * PE # Update Q-values
      }  
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M03: RL-DDM with separate bias terms for Win and Avoid cues: #### 

modSim_M03 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial beta:
      beta <- par[iSub, 3] * data$valence[iTrialSub] + par[iSub, 4] * (1 - data$valence[iTrialSub])
      
      ## Compute trial-by-trial delta:
      delta = par[iSub, 5] + (Q[s, 1] - Q[s, 2]) * par[iSub, 6]
      
      ## Simulate response and RT using rwiener:
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = beta, delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      } 
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M04: RL-DDM with separate drift intercept terms for Win and Avoid cues: #### 

modSim_M04 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M05: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold term for high/low stakes: #### 

modSim_M05 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      alpha <- par[iSub, 1] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = alpha, tau = par[iSub, 2], beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M06: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time for high/low stakes: #### 

modSim_M06 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      tau <- par[iSub, 2] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = tau, beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M07: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate bias term for high/low stakes: #### 

modSim_M07 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      beta <- par[iSub, 3] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = beta, delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M08: RL-DDM with separate drift intercept terms for Win and Avoid cues and drift bonus for high stakes: #### 

modSim_M08 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub]) + 
        par[iSub, 8] * data$stakes[iTrialSub]
      
      ## Simulate response and RT using rwiener:
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = par[iSub, 2], beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M09: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold and non-decision time for high/low stakes: #### 

modSim_M09 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      alpha <- par[iSub, 1] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      tau <- par[iSub, 2] * (1 - data$stakes[iTrialSub]) + par[iSub, 9] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = alpha, tau = tau, beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M10: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold and drift rate bonus for high/low stakes: #### 

modSim_M10 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub]) + 
        par[iSub, 9] * data$stakes[iTrialSub]
      
      ## Simulate response and RT using rwiener:
      alpha <- par[iSub, 1] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = alpha, tau = par[iSub, 2], beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M11: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time and drift bonus for high/low stakes: #### 

modSim_M11 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub]) + 
        par[iSub, 9] * data$stakes[iTrialSub]
      
      ## Simulate response and RT using rwiener:
      tau <- par[iSub, 2] * (1 - data$stakes[iTrialSub]) + par[iSub, 8] * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = tau, beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}

# ============================================================================================================================= #
#### modSim_M12: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time for low, high/congruent, and high/incongruent stakes: #### 

modSim_M12 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables subject, trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nSub x nPar matrix with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nSub <- length(unique(data$subject))
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  for (iSub in 1:nSub){
    
    subIdx <- which(data$subject == iSub)
    # table(data$valence[subIdx], data$stimuli[subIdx])
    valPerStim <- tapply(data$valence[subIdx], data$stimuli[subIdx], mean) - 0.5
    Q <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
    
    for (iTrial in 1:nTrial){
      
      ## Extract trial data:
      iTrialSub <- (iSub-1) * nTrial + iTrial
      s <- data$stimuli[iTrialSub]
      c <- data$resp[iTrialSub]
      r <- data$outcome[iTrialSub]
      
      data$QGo[iTrialSub] <- Q[s, 1]
      data$QNoGo[iTrialSub] <- Q[s, 2]
      
      ## Compute trial-by-trial delta:
      delta <- (Q[s, 1] - Q[s, 2]) * par[iSub, 4] + 
        par[iSub, 5] * data$valence[iTrialSub] + 
        par[iSub, 6] * (1 - data$valence[iTrialSub])
      
      ## Simulate response and RT using rwiener:
      tau <- par[iSub, 2] * (1 - data$stakes[iTrialSub]) + 
        (par[iSub, 8] * data$congruency[iTrialSub] + par[iSub, 9] * (1 - data$congruency[iTrialSub])) * data$stakes[iTrialSub]
      dat <- rwiener(n = 1, alpha = par[iSub, 1], tau = tau, beta = par[iSub, 3], delta = delta)
      data$sim_resp[iTrialSub] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
      data$sim_rt[iTrialSub] <- dat$q
      c <- data$sim_resp[iTrialSub] # save, crop too slow resps for updating
      
      ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
      if(c == 1 & dat$q > 1.3){c <- 0}
      
      ## Simulate outcome:
      r <- NA # default
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] == c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- -1}
      
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 1) {r <- 0}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 1 & data$valence[iTrialSub] == 0) {r <- -1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 1) {r <- 1}
      if(data$reqAction[iTrialSub] != c & data$validity[iTrialSub] == 0 & data$valence[iTrialSub] == 0) {r <- 0}
      
      data$sim_outcome[iTrialSub] <- r
      
      ## Update Q-values:
      if (!is.na(r)) { # update only if outcome not NA 
        PE <- r - Q[s, 2-c] # Compute prediction error
        Q[s, 2-c] <- Q[s, 2-c] + par[iSub, 7] * PE # Update Q-values
      }
      
    } # end iTrial
  } # end iSub
  return(data)
}



# END OF FILE.