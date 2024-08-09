# RLDDM4stan_sepParam_flat.R
# ============================================================================================================================= #
# Specification model strings for rstan.
# Drift rate intercept AND slope.
# Effect of stakes modeled as separate parameter (not difference parameter) with same initialization/transformation/prior as original parameter.
# adapted by Johannes Algermissen, 2023.

# M01: DDM with one drift rate (chance performance with bias).
# M02: RL-DDM (needs cue, response, outcome).
# M03: RL-DDM with separate starting points.
# M04: RL-DDM with separate drift rates.

# consider fixing threshold to 1
# consider using stakes differences or dwell time difference/ ratio
# consider non-linear drift rate (see Fontanesi) 

# Uses wiener_lpdf: source code under https://mc-stan.org/math/d5/d2f/wiener__lpdf_8hpp_source.html

# ============================================================================================================================= #
#### mstr01: DDM (overall pGo): ####
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), delta (drift rate).
# With one unconstrained drift rate (can be both positive or negative)
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierDDM.stan
# ============================================================================================================================= #

# DO NOT flip beta and delta for lower bound; handle with wiener_gonogo_lpdf

mstr01 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nTrial;                 // total number of trials
  
  int<lower=-1,upper=1> resp[nTrial];  // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                   // RTs (in sec.)

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaInt;                   // delta: drift rate

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real deltaInt;                      // delta: drift rate (both positive and negative)

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);
  tau = log1p_exp(nd_tau);
  beta = inv_logit(nd_beta);          // starting point bias bound to range 0-1
  deltaInt = nd_deltaInt;             // do not transform because must be allowed to positive or negative

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.9, 1);                // alpha: boundary separation/threshold
  // nd_tau ~ normal(-2.0, 1);                 // tau: nondecision time
  // nd_beta ~ normal(-2.7, 1);                // beta: starting point bias
  // nd_deltaInt ~ normal(3.6, 1);             // delta: drift rate
  
  // Informative priors:
  nd_alpha ~ normal(1.87, 0.31);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-2.00, 0.60);           // tau: nondecision time
  nd_beta ~ normal(-2.74, 0.35);          // beta: starting point bias
  nd_deltaInt ~ normal(3.62, 0.41);       // delta: drift rate

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaInt);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaInt);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr02: RL-DDM: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr02 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaInt;                   // deltaInt: drift rate intercept (Go bias)
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_epsilon;                    // epsilon: learning rate

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real deltaInt;                      // deltaInt: drift rate intercept (Go bias)
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaInt = nd_deltaInt;                 // Go bias, no transform
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaInt + deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);                // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);                 // tau: nondecision time
  // nd_beta ~ normal(-1.1, 1);                // beta: starting point bias
  // nd_deltaInt ~ normal(1.3, 1);             // deltaInt: drift rate intercept (Go bias)
  // nd_deltaSlope ~ normal(1.7, 1);           // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_epsilon ~ normal(0, 1);                // epsilon: learning rate
  
  // Informative priors:
  nd_alpha ~ normal(1.08, 0.21);            // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.33, 0.40);             // tau: nondecision time
  nd_beta ~ normal(-1.10, 0.37);            // beta: starting point bias
  nd_deltaInt ~ normal(1.28, 0.45);         // deltaInt: drift rate intercept (Go bias)
  nd_deltaSlope ~ normal(1.70, 0.81);       // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_epsilon ~ normal(-0.05, 0.69);         // epsilon: learning rate

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr03: RL-DDM with separate starting points: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# betaWin/ betaAvoid (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr03 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_betaWin;                    // betaWin: starting point bias for Win cues
  real nd_betaAvoid;                  // betaAvoid: starting point bias for Avoid cues
  real nd_deltaInt;                   // deltaInt: drift rate intercept (Go bias)
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_epsilon;                    // epsilon: learning rate

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> betaWin;    // betaWin: starting point bias for Win cues
  real <lower=0, upper=1> betaAvoid;  // betaAvoid: starting point bias for Avoid cues
  real deltaInt;                      // deltaInt: drift rate intercept (Go bias)
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  
  vector[nTrial] betaTrial;           // betaTrial: trial-by-trial starting point bias (depending on cue valence)
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  betaWin = inv_logit(nd_betaWin);        // starting point bias bound to range 0-1
  betaAvoid = inv_logit(nd_betaAvoid);    // starting point bias bound to range 0-1
  deltaInt = nd_deltaInt;                 // Go bias, no transform
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials

    // Compute trial-by-trial starting point bias:
    betaTrial[iTrial] = betaWin * valence[iTrial] + betaAvoid * (1 - valence[iTrial]);

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaInt + deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2]);

    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);              // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.4, 1);               // tau: nondecision time
  // nd_betaWin ~ normal(-1.0, 1);           // betaWin: starting point bias for Win cues
  // nd_betaAvoid~ normal(-1.5, 1);          // betaAvoid: starting point bias for Win cues
  // nd_deltaInt ~ normal(1.5, 1);           // deltaInt: drift rate intercept (Go bias)
  // nd_deltaSlope ~ normal(1.6, 1);         // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_epsilon ~ normal(0, 1);              // epsilon: learning rate
  
  // Informative priors:
  nd_alpha ~ normal(1.14, 0.23);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.36, 0.42);           // tau: nondecision time
  nd_betaWin ~ normal(-0.99, 0.32);       // betaWin: starting point bias for Win cues
  nd_betaAvoid ~ normal(-1.53, 0.43);     // betaAvoid: starting point bias for Avoid cues
  nd_deltaInt ~ normal(1.48, 0.59);       // deltaInt: drift rate intercept (Go bias)
  nd_deltaSlope ~ normal(1.56, 0.96);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_epsilon ~ normal(0.00, 0.65);        // epsilon: learning rate

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, betaTrial[iTrial], deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, betaTrial[iTrial], deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr04: RL-DDM with separate reward and punishment modulation on drift: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaWin/ deltaAvoid (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr04 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.1, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.6, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  
  // Informative priors:
  nd_alpha ~ normal(1.10, 0.19);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.33, 0.39);           // tau: nondecision time
  nd_beta ~ normal(-1.07, 0.31);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.51, 0.91);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.64, 0.60);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.92, 0.33);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.13, 0.58);        // epsilon: learning rate

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr05: RL-DDM with separate reward and punishment modulation on drift and and separate threshold for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate threshold for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr05 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate threshold for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate threshold for high stakes
  
  vector[nTrial] alphaTrial;          // alphaTrial: trial-by-trial decision threshold
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // decision threshold, must be positive

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial alpha:
    alphaTrial[iTrial] = pi*stakes[iTrial] + alpha*(1-stakes[iTrial]);
    
    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.1, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.7, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(1.1, 1);               // pi: separate threshold for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.05, 0.18);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.32, 0.38);           // tau: nondecision time
  nd_beta ~ normal(-1.07, 0.31);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.52, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.65, 0.61);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.93, 0.33);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.13, 0.59);        // epsilon: learning rate
  nd_pi ~ normal(1.14, 0.18);             // pi: separate threshold for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr06: RL-DDM with separate reward and punishment modulation on drift and separate non-decision time for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate non-decision time for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr06 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate nondecision time for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate nondecision time for high stakes
  
  vector[nTrial] tauTrial;            // tauTrial: trial-by-trial non-decision time
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // non-decision time, must be positive

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial tau:
    tauTrial[iTrial] = pi*stakes[iTrial] + tau*(1-stakes[iTrial]);

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.0, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.6, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(-1.2, 1);              // pi: nondecision time for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.05, 0.17);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.30, 0.36);           // tau: nondecision time
  nd_beta ~ normal(-1.00, 0.29);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.54, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.58, 0.61);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.85, 0.31);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.12, 0.60);        // epsilon: learning rate
  nd_pi ~ normal(-1.24, 0.36);            // pi: nondecision time for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr07: RL-DDM with separate reward and punishment modulation on drift and separate bias for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (effect of stakes on bias).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr07 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate bias for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0, upper=1> pi;         // pi: separate bias for high stakes

  vector[nTrial] betaTrial;           // betaTrial: trial-by-trial starting point bias
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = inv_logit(nd_pi);                  // starting point bias bound to range 0-1

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial beta:
    betaTrial[iTrial] = pi*stakes[iTrial] + beta*(1-stakes[iTrial]); 

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.0, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.6, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(-1.1, 1);              // pi: starting point bias for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.09, 0.18);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.32, 0.38);           // tau: nondecision time
  nd_beta ~ normal(-0.99, 0.23);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.51, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.61, 0.64);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.88, 0.33);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.13, 0.59);        // epsilon: learning rate
  nd_pi ~ normal(-1.08, 0.21);            // pi: starting point bias for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, betaTrial[iTrial], deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, betaTrial[iTrial], deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr08: RL-DDM with separate reward and punishment modulation on drift and drift bonus for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (drift bonus for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr08 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: drift bonus for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real pi;                            // pi: drift bonus for high stakes

  vector[nTrial] deltaTrial;          // deltaTrial: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = nd_pi;                             // drift bonus for high stakes

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial])
    + pi * stakes[iTrial];
 
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.1, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.7, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(1.0, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(-0.1, 1);              // pi: drift bonus for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.10, 0.19);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.33, 0.39);           // tau: nondecision time
  nd_beta ~ normal(-1.08, 0.31);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.51, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.70, 0.61);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.98, 0.33);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.13, 0.59);        // epsilon: learning rate
  nd_pi ~ normal(-0.10, 0.04);            // pi: drift bonus for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr09: RL-DDM with separate reward and punishment modulation on drift and separate threshold and non-decision time for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate threshold for high stakes), theta (separate non-decision time for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr09 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate threshold for high stakes
  real nd_theta;                      // theta: separate nondecision time for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate threshold for high stakes
  real <lower=0> theta;               // theta: separate nondecision time for high stakes
  
  vector[nTrial] alphaTrial;          // alphaTrial: trial-by-trial decision threshold
  vector[nTrial] tauTrial;            // tauTrial: trial-by-trial non-decision time
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // decision threshold, must be positive
  theta = log1p_exp(nd_theta);            // non-decision time, must be positive

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial alpha:
    alphaTrial[iTrial] = pi*stakes[iTrial] + alpha*(1-stakes[iTrial]);

    // Compute trial-by-trial tau:
    tauTrial[iTrial] = theta*stakes[iTrial] + tau*(1-stakes[iTrial]);
    
    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.0, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.0, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.6, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.8, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(1.1, 1);               // pi: separate threshold for high stakes
  // nd_theta ~ normal(-1.3, 1);           // theta: nondecision time
  
  // Informative priors:
  nd_alpha ~ normal(1.01, 0.19);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.27, 0.37);           // tau: nondecision time
  nd_beta ~ normal(-0.98, 0.29);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.54, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.58, 0.60);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.84, 0.31);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.12, 0.60);        // epsilon: learning rate
  nd_pi ~ normal(1.08, 0.19);             // pi: separate threshold for high stakes
  nd_theta ~ normal(-1.25, 0.36);         // theta: nondecision time

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr10: RL-DDM with separate reward and punishment modulation on drift and separate threshold and drift bonus for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate threshold for high stakes), theta (drift bonus for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr10 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate threshold for high stakes
  real nd_theta;                      // theta: drift bonus for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate threshold for high stakes
  real theta;                         // theta: drift bonus for high stakes
  
  vector[nTrial] alphaTrial;          // alphaTrial: trial-by-trial decision threshold
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // decision threshold, must be positive
  theta = nd_theta;                       // drift bonus, no transform

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial alpha:
    alphaTrial[iTrial] = pi*stakes[iTrial] + alpha*(1-stakes[iTrial]);

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial])
    + theta * stakes[iTrial];
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.1, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.7, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(1.1, 1);               // pi: separate threshold for high stakes
  // nd_theta ~ normal(0.0, 1);            // theta: drift bonus for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.05, 0.18);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.32, 0.38);           // tau: nondecision time
  nd_beta ~ normal(-1.07, 0.31);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.52, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.67, 0.61);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.95, 0.33);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.13, 0.59);        // epsilon: learning rate
  nd_pi ~ normal(1.14, 0.18);             // pi: separate threshold for high stakes
  nd_theta ~ normal(-0.03, 0.04);         // theta: drift bonus for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alphaTrial[iTrial], tau, beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr11: RL-DDM with separate reward and punishment modulation on drift and separate non-decision time and drift bonus for high stakes: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate non-decision time for high stakes), theta (drift bonus for high stakes).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr11 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate nondecision time for high stakes
  real nd_theta;                      // theta: drift bonus for high stakes

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate nondecision time for high stakes
  real theta;                         // theta: drift bonus for high stakes
  
  vector[nTrial] tauTrial;            // tauTrial: trial-by-trial nondecision time
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // non-decision time, must be positive
  theta = nd_theta;                       // drift bonus, no transform

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial alpha:
    tauTrial[iTrial] = pi*stakes[iTrial] + tau*(1-stakes[iTrial]);

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial])
    + theta * stakes[iTrial];
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-1.0, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.5, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.6, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.9, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(-1.2, 1);              // pi: separate nondecision time for high stakes
  // nd_theta ~ normal(-0.1, 1);           // theta: drift bonus for high stakes
  
  // Informative priors:
  nd_alpha ~ normal(1.05, 0.17);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.29, 0.36);           // tau: nondecision time
  nd_beta ~ normal(-1.00, 0.29);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.54, 0.92);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.61, 0.60);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.88, 0.31);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.12, 0.60);        // epsilon: learning rate
  nd_pi ~ normal(-1.25, 0.36);            // pi: separate nondecision time for high stakes
  nd_theta ~ normal(-0.05, 0.05);         // theta: drift bonus for high stakes

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# ============================================================================================================================= #
#### mstr12: RL-DDM with separate reward and punishment modulation on drift and separate non-decision time for high stakes separately for congruent/incongruent cues: #### 
# Contains parameters alpha (boundary separation/threshold), tau (non-decision time),
# beta (starting point bias), deltaInt (Go bias), deltaSlope (scaling of Q-value difference), 
# epsilon (learning rate), pi (separate non-decision time for high stakes for congruent cues), theta (separate non-decision time for high stakes for incongruent cues).
# after https://github.com/laurafontanesi/rlssm/blob/main/rlssm/stan_models/hierRLDDM.stan
# ============================================================================================================================= #

mstr12 <- "
functions {
  /* returns Wiener First Passage Time Distribution if resp == 1 (Go), i.e., when RT is available
   * otherwise integrates over RTs (NoGo) based on solution in Blurton, Kesselmeier, & Gondan (2012, Journal of Mathematical Psychology) 
   * see also discussion under https://groups.google.com/g/hddm-users/c/3hlU6RnNcrw ,
   * and implementation in function prob_ub in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/pdf.pxi ,
   * which is called in function wiener_like_multi in HDDM: https://github.com/hddm-devs/hddm/blob/master/src/wfpt.pyx .
   * Args:
   *   RT   : vector or reaction times (in seconds)
   *   resp : vector or responses, either 1 (Go) or 0 (NoGo)
   *   alpha: Boundary separation
   *   tau  : Nondecision time
   *   beta : A-priori bias
   *   delta: Drift rate
   * Returns:
   *   log probability density for RT/resp given parameters
   */
  real wiener_gonogo_lpdf(real RT, real resp, real alpha, real tau, real beta, real delta) { 
    if (resp == 1) { // Go: just use wiener_lpdf
      return wiener_lpdf(RT | alpha, tau, beta, delta);
    } else if (delta == 0) { // NoGo and drift == 0
      return log1m(beta); // just return starting point bias (1 minus because lower boundary)
    } else { // NoGo and drift != 0
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // numerical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of stimuli
  int<lower=0> nResp;                     // total number of responses
  int<lower=0> nTrial;                    // total number of trials
  
  int<lower=1> stimuli[nTrial];           // cue identifier
  int<lower=-1,upper=1> resp[nTrial];     // choices (1=upper bound, 0=lower bound)
  vector[nTrial] rt;                      // RTs (in sec.)
  
  int<lower=-1, upper=1> outcome[nTrial]; // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nTrial];  // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nTrial];   // stakes (1=high, 0=low)
  int<lower=0, upper=1> congruency[nTrial];// congruency (1=congruent, 0=incongruent)

  matrix[nStim,nResp] Qi;                 // Initial action values for all cues: for each cue, for each response

}
parameters {    

  // Define parameters:
  real nd_alpha;                      // alpha: boundary separation (threshold)
  real nd_tau;                        // tau: nondecision time
  real nd_beta;                       // beta: starting point bias
  real nd_deltaSlope;                 // deltaSlope: drift rate slope (scaling Q-value difference)
  real nd_deltaWin;                   // deltaWin: drift rate offset for Win cues
  real nd_deltaAvoid;                 // deltaAvoid: drift rate offset for Avoid cues
  real nd_epsilon;                    // epsilon: learning rate
  real nd_pi;                         // pi: separate nondecision time for high stakes for congruent cues
  real nd_theta;                      // theta: separate nondecision time for high stakes for incongruent cues

} // end parameters
transformed parameters {

  // Define transformed parameters:
  real <lower=0> alpha;               // alpha: boundary separation (threshold)
  real <lower=0> tau;                 // tau: nondecision time
  real <lower=0, upper=1> beta;       // beta: starting point bias
  real <lower=0> deltaSlope;          // deltaSlope: drift rate (scaling Q-value difference)
  real deltaWin;                      // deltaWin: drift rate offset for Win cues
  real deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> epsilon;    // epsilon: learning rate
  real <lower=0> pi;                  // pi: separate nondecision time for high stakes for congruent cues
  real <lower=0> theta;               // theta: separate nondecision time for high stakes for incongruent cues
  
  vector[nTrial] tauTrial;            // tauTrial: trial-by-trial non-decision time
  vector[nTrial] deltaTrial;          // deltaTrial: trial-by-trial drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;              // Current action values for all cues: for each cue, for each response
  real PE;                            // Prediction error on this trial

  // Transform parameters:
  alpha = log1p_exp(nd_alpha);            // decision threshold, must be positive
  tau = log1p_exp(nd_tau);                // non-decision time, must be positive
  beta = inv_logit(nd_beta);              // starting point bias bound to range 0-1
  deltaSlope = log1p_exp(nd_deltaSlope);  // transform to be positive (do not invert QGo-QNoGo)
  deltaWin = nd_deltaWin;                 // Go bias, no transform
  deltaAvoid = nd_deltaAvoid;             // Go bias, no transform
  epsilon = inv_logit(nd_epsilon);        // learning rate bound to range 0-1 
  pi = log1p_exp(nd_pi);                  // non-decision time, must be positive
  theta = log1p_exp(nd_theta);            // non-decision time, must be positive

  // Initialize variables:
  Q = Qi;                                 // initialize action values (stored in simulated data; copy over)

  for(iTrial in 1:nTrial){                // Loop over trials
    
    // Compute trial-by-trial tau:
    tauTrial[iTrial] = (pi * congruency[iTrial] + theta * (1-congruency[iTrial])) * stakes[iTrial] + tau * (1-stakes[iTrial]);

    // Compute trial-by-trial drift rate:
    deltaTrial[iTrial] = deltaSlope * (Q[stimuli[iTrial], 1] - Q[stimuli[iTrial], 2])
    + deltaWin * valence[iTrial]
    + deltaAvoid * (1 - valence[iTrial]);
      
    // Update Q-values:
    PE = outcome[iTrial] - Q[stimuli[iTrial], 2-resp[iTrial]];                            // Compute prediction error
    Q[stimuli[iTrial], 2-resp[iTrial]] = Q[stimuli[iTrial], 2-resp[iTrial]] + epsilon*PE; // Update Q-values
      
  } // end iTrial

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Uninformative priors:
  // nd_alpha ~ normal(1.1, 1);            // alpha: boundary separation/threshold
  // nd_tau ~ normal(-1.3, 1);             // tau: nondecision time
  // nd_beta ~ normal(-0.9, 1);            // beta: starting point bias
  // nd_deltaSlope ~ normal(1.6, 1);       // deltaSlope: drift rate slope (scaling Q-value difference)
  // nd_deltaWin ~ normal(1.5, 1);         // deltaWin: drift rate offset for Win cues
  // nd_deltaAvoid ~ normal(0.8, 1);       // deltaAvoid: drift rate offset for Avoid cues
  // nd_epsilon ~ normal(0.1, 1);          // epsilon: learning rate
  // nd_pi ~ normal(-1.2, 1);              // pi: nondecision time for high stakes for congruent cues
  // nd_theta ~ normal(-1.2, 1);           // theta: nondecision time for high stakes for incongruent cues
  
  // Informative priors:
  nd_alpha ~ normal(1.01, 0.16);          // alpha: boundary separation/threshold
  nd_tau ~ normal(-1.26, 0.35);           // tau: nondecision time
  nd_beta ~ normal(-0.92, 0.29);          // beta: starting point bias
  nd_deltaSlope ~ normal(1.58, 0.90);     // deltaSlope: drift rate slope (scaling Q-value difference)
  nd_deltaWin ~ normal(1.52, 0.61);       // deltaWin: drift rate offset for Win cues
  nd_deltaAvoid ~ normal(0.79, 0.32);     // deltaAvoid: drift rate offset for Avoid cues
  nd_epsilon ~ normal(0.11, 0.64);        // epsilon: learning rate
  nd_pi ~ normal(-1.17, 0.30);            // pi: nondecision time for high stakes for congruent cues
  nd_theta ~ normal(-1.18, 0.35);         // theta: nondecision time for high stakes for incongruent cues

  // Loop over data, increment log-likelihood:
  for (iTrial in 1:nTrial) {
      target += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iTrial in 1:nTrial) {
      log_lik += wiener_gonogo_lpdf( rt[iTrial] | resp[iTrial], alpha, tauTrial[iTrial], beta, deltaTrial[iTrial]);
  } // end iTrial 
  
} // end generated quantities
"

# END OF FILE.
