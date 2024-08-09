#### Model simulations single subject: ####
# simulate response and outcomes anew
# only for a single subject

# ============================================================================================================================= #
#### modSim_M01: DDM (chance performance with bias) ####

modSim_M01 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA # will stay NA
  data$QNoGo <- NA # will stay NA
  data$resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Simulate response and RT using rwiener:
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = par[3], delta = par[4])
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save, crop too slow responses for updating
    
    # ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
  } # end iTrial
  
  ## Save data as list:
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
  
}

# ============================================================================================================================= #
#### modSim_02: RL-DDM (learning): #### 

modSim_M02 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
    
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = par[4] + (Q[s, 1] - Q[s, 2]) * par[5]
    
    ## Simulate response and RT using rwiener:
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[6] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M03: RL-DDM with separate bias terms for Win and Avoid cues: #### 

modSim_M03 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial beta:
    beta <- par[3] * data$valence[iTrial] + par[4] * (1 - data$valence[iTrial])
    
    ## Compute trial-by-trial delta:
    delta = par[5] + (Q[s, 1] - Q[s, 2]) * par[6]
    
    ## Simulate response and RT using rwiener:
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = beta, delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}
# ============================================================================================================================= #
#### modSim_M04: RL-DDM with separate drift intercept terms for Win and Avoid cues: #### 

modSim_M04 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M05: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold term for high/low stakes: #### 

modSim_M05 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    alpha <- par[1] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = alpha, tau = par[2], beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M06: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time for high/low stakes: #### 

modSim_M06 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    tau <- par[2] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = par[1], tau = tau, beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M07: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate bias term for high/low stakes: #### 

modSim_M07 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    beta <- par[3] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = beta, delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M08: RL-DDM with separate drift intercept terms for Win and Avoid cues and drift bonus for high stakes: #### 

modSim_M08 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial]) + 
      par[8] * data$stakes[iTrial]
    
    ## Simulate response and RT using rwiener:
    dat <- rwiener(n = 1, alpha = par[1], tau = par[2], beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M09: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold and non-decision time for high/low stakes: #### 

modSim_M09 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    alpha <- par[1] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    tau <- par[2] * (1 - data$stakes[iTrial]) + par[9] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = alpha, tau = tau, beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M10: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate threshold and drift rate bonus for high/low stakes: #### 

modSim_M10 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial]) +
      par[9] * data$stakes[iTrial]
    
    ## Simulate response and RT using rwiener:
    alpha <- par[1] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = alpha, tau = par[2], beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M11: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time and drift bonus for high/low stakes: #### 

modSim_M11 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial]) + 
      par[9] * data$stakes[iTrial]
    
    ## Simulate response and RT using rwiener:
    tau <- par[2] * (1 - data$stakes[iTrial]) + par[8] * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = par[1], tau = tau, beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# ============================================================================================================================= #
#### modSim_M12: RL-DDM with separate drift intercept terms for Win and Avoid cues and separate non-decision time for low, high/congruent, and high/incongruent stakes: #### 

modSim_M12 <- function (data, par){
  #' Simulate responses and RTs based on input data and fitted parameters
  #' data   - data frame with variables trialnr, stimuli, reqAction, valence, stimulation, validity, resp, rt, outcome
  #' par    - nPar vector with subject-level parameters
  #' output: input data with columns sim_resp and sim_rt and sim_outcome added
  
  require(RWiener)
  
  nTrial <- length(unique(data$trialnr))
  nStim <- length(unique(data$stimuli))
  nResp <- length(unique(data$resp))
  
  ## Initialize empty variables:
  data$QGo <- NA
  data$QNoGo <- NA
  data$sim_resp <- NA
  data$sim_rt <- NA
  
  valPerStim <- tapply(data$valence, data$stimuli, mean) - 0.5
  Qi <- matrix(valPerStim, nrow = nStim, ncol = nResp, byrow = F) # initialize Q-values
  Q <- Qi
  
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    ## Extract trial data:
    s <- data$stimuli[iTrial]
    c <- data$resp[iTrial]
    r <- data$outcome[iTrial]
    
    data$QGo[iTrial] <- Q[s, 1]
    data$QNoGo[iTrial] <- Q[s, 2]
    
    ## Compute trial-by-trial delta:
    delta = (Q[s, 1] - Q[s, 2]) * par[4] + 
      par[5] * data$valence[iTrial] + 
      par[6] * (1 - data$valence[iTrial])
    
    ## Simulate response and RT using rwiener:
    tau <- par[2] * (1 - data$stakes[iTrial]) + 
      (par[8] * data$congruency[iTrial] + par[9] * (1 - data$congruency[iTrial])) * data$stakes[iTrial]
    dat <- rwiener(n = 1, alpha = par[1], tau = tau, beta = par[3], delta = delta)
    data$sim_resp[iTrial] <- ifelse(dat$resp == "upper", 1, 0) # save even late resps, crop later for updating
    data$sim_rt[iTrial] <- dat$q
    c <- data$sim_resp[iTrial] # save; crop too slow resps for updating
    
    ## Set c to NoGo if RT after deadline (but still save resp and rt as Go response above):
    if(c == 1 & dat$q > 1.3){c <- 0}
    
    ## Simulate outcome:
    r <- NA # default
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] == c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- -1}
    
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 1) {r <- 0}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 1 & data$valence[iTrial] == 0) {r <- -1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 1) {r <- 1}
    if(data$reqAction[iTrial] != c & data$validity[iTrial] == 0 & data$valence[iTrial] == 0) {r <- 0}
    
    data$sim_outcome[iTrial] <- r
    
    ## Update Q-values:
    if (!is.na(r)) { # update only if outcome not NA 
      PE <- r - Q[s, 2-c] # Compute prediction error
      Q[s, 2-c] <- Q[s, 2-c] + par[7] * PE # Update Q-values
    }  
    
  } # end iTrial
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, trialnr = data$trialnr,
                stimuli = data$stimuli, reqAction = data$reqAction, valence = data$valence, congruency = data$congruency, stakes = data$stakes,
                resp = data$sim_resp, rt = data$sim_rt, validity = data$validity, outcome = data$sim_outcome, 
                RTbound = 0, minRT = min(data$sim_rt), Qi = Qi)
  return(dList)
}

# END OF FILE.