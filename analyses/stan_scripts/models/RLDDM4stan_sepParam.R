# RLDDM4stan_sepParam.R
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nTrial;              // number of trials
  int<lower=0> nSub;                // number of subject
  int<lower=0> nData;               // total number of data points
  
  int<lower=0,upper=1> resp[nData]; // choices (1=upper bound, 0=lower bound)
  vector[nData] rt;                 // RTs (in sec.)

  real RTbound;                     // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                       // minimum RT of the observed data (upper limit for tau)

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaInt;                             // delta: drift rate

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaInt;          // delta: drift rate

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaInt[nSub];                        // delta: drift rate

} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaInt;                      // delta: drift rate

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real sub_deltaInt[nSub];                      // delta: drift rate (both positive and negative)

  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaInt = mu_deltaInt; //  don't transform, but keep

  for (iSub in 1:nSub){
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);      // starting point bias bound to range 0-1
    sub_deltaInt[iSub] = mu_deltaInt + z_deltaInt[iSub] * sd_deltaInt;  // don't transform because must be allowed to positive or negative

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 
    delta[firstTrial:lastTrial] = rep_vector(sub_deltaInt[iSub],nTrial); 
  }
  
} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaInt ~ normal(0, 3);       // delta: drift rate (intercept)

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaInt ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaInt ~ normal(0, 1);

  // Loop over subject, loop over trials, increment log-likelihood:
  for (iData in 1:nData) {
      target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData 
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
      log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData 
  
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaInt;                             // deltaInt: drift rate intercept (Go bias)
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_epsilon;                              // epsilon: learning rate

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaInt;          // deltaInt: drift rate intercept (Go bias)
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaInt[nSub];                        // deltaInt: drift rate intercept (Go bias)
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_epsilon[nSub];                         // epsilon: learning rate
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaInt;                      // deltaInt: drift rate intercept (Go bias)
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_epsilon;                       // epsilon: learning rate

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real sub_deltaInt[nSub];                      // deltaInt: drift rate intercept (Go bias)
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaInt = mu_deltaInt;
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_epsilon = inv_logit(mu_epsilon);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                          // starting point bias bound to range 0-1
    sub_deltaInt[iSub] = mu_deltaInt + z_deltaInt[iSub] * sd_deltaInt;                      // Go bias, no transform
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // transform to be positive (don't invert QGo-QNoGo)
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);              // learning rate bound to range 0-1

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];          // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];          // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaInt[iSub] + sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]);
        
      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaInt ~ normal(0, 3);       // deltaInt: drift rate intercept (Go bias)
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaInt ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaInt ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)

  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_betaWin;                              // betaWin: starting point bias for Win cues
  real mu_betaAvoid;                            // betaAvoid: starting point bias for Avoid cues
  real mu_deltaInt;                             // deltaInt: drift rate intercept (Go bias)
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_epsilon;                              // epsilon: learning rate

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_betaWin;           // betaWin: starting point bias for Win cues
  real<lower=0, upper=20> sd_betaAvoid;         // betaAvoid: starting point bias for Avoid cues
  real<lower=0, upper=20> sd_deltaInt;          // deltaInt: drift rate intercept (Go bias)
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_betaWin[nSub];                         // betaWin: starting point bias for Win cues
  real z_betaAvoid[nSub];                       // betaAvoid: starting point bias for Avoid cues
  real z_deltaInt[nSub];                        // deltaInt: drift rate intercept (Go bias)
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_epsilon[nSub];                         // epsilon: learning rate
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_betaWin;                       // betaWin: starting point bias for Win cues
  real transf_mu_betaAvoid;                     // betaAvoid: starting point bias for Avoid cues
  real transf_mu_deltaInt;                      // deltaInt: drift rate intercept (Go bias)
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_epsilon;                       // epsilon: learning rate

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_betaWin[nSub];    // betaWin: starting point bias for Win cues
  real <lower=0, upper=1> sub_betaAvoid[nSub];  // betaAvoid: starting point bias for Avoid cues
  real sub_deltaInt[nSub];                      // deltaInt: drift rate intercept (Go bias)
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_betaWin = inv_logit(mu_betaWin);
  transf_mu_betaAvoid = inv_logit(mu_betaAvoid);
  transf_mu_deltaInt = mu_deltaInt;
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_epsilon = inv_logit(mu_epsilon);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);
    sub_betaWin[iSub] = inv_logit(mu_betaWin + z_betaWin[iSub] * sd_betaWin);              // starting point bias for Win cues bound to range 0-1
    sub_betaAvoid[iSub] = inv_logit(mu_betaAvoid + z_betaAvoid[iSub] * sd_betaAvoid);      // starting point bias for Avoid cues bound to range 0-1
    sub_deltaInt[iSub] = mu_deltaInt + z_deltaInt[iSub] * sd_deltaInt;                      // Go bias, no transform
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // transform to be positive (don't invert QGo-QNoGo)
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);              // learning rate bound to range 0-1

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial);

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial starting point:
      beta[iTrialSub] = sub_betaWin[iSub] * valence[iTrialSub] + sub_betaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaInt[iSub] + sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]);
        
      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_betaWin ~ normal(0, 1);        // betaWin: starting point bias for Win cues
  mu_betaAvoid ~ normal(0, 1);      // betaAvoid: starting point bias for Avoid cues
  mu_deltaInt ~ normal(0, 3);       // deltaInt: drift rate intercept (Go bias)
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_betaWin ~ normal(0, 1);
	sd_betaAvoid ~ normal(0, 1);
	sd_deltaInt ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_betaWin ~ normal(0, 1);
	z_betaAvoid ~ normal(0, 1);
	z_deltaInt ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                          // starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);              // learning rate bound to range 0-1

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate threshold for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate threshold for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate threshold for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate threshold for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate threshold for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: transform to be positive

    // Repeat for each trial of this subject:
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub], nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub], nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial alpha:
      alpha[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_alpha[iSub]*(1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate threshold for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate non-decision time for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate non-decision time for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate non-decision time for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate non-decision time for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate non-decision time for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: transform to be positive

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub], nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub], nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial tau:
      tau[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_tau[iSub]*(1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate non-decision time for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes(1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate bias for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate bias for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate bias for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate bias for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0, upper=1> sub_pi[nSub];         // pi: separate bias for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = inv_logit(mu_pi);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = inv_logit(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: separate bias for high stakes

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial);

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial beta:
      beta[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_beta[iSub]*(1-stakes[iTrialSub]); 

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate bias for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes(1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: drift bonus for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: drift bonus for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: drift bonus for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: drift bonus for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real sub_pi[nSub];                            // pi: drift bonus for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = mu_pi;
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = mu_pi + z_pi[iSub] * sd_pi;                                              // pi: drift bonus for high stakes

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub], nTrial);
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub], nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub], nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub], 1] - Q[stimuli[iTrialSub], 2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]) 
      + sub_pi[iSub] * stakes[iTrialSub];

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub], 2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub], 2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: drift bonus for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate threshold for high stakes
  real mu_theta;                                // theta: separate non-decision time for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate threshold for high stakes
  real<lower=0, upper=20> sd_theta;             // theta: separate non-decision time for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate threshold for high stakes
  real z_theta[nSub];                           // theta: separate non-decision time for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate threshold for high stakes
  real transf_mu_theta;                         // theta: separate non-decision time for high stakes
  
  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate threshold for high stakes
  real <lower=0> sub_theta[nSub];               // theta: separate non-decision time for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  transf_mu_theta = log1p_exp(mu_theta);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: separate threshold for high stakes
    sub_theta[iSub] = log1p_exp(mu_theta + z_theta[iSub] * sd_theta);                       // theta: separate non-decision time for high stakes

    // Repeat for each trial of this subject:
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial alpha and tau:
      alpha[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_alpha[iSub]*(1-stakes[iTrialSub]);
      tau[iTrialSub] = sub_theta[iSub]*stakes[iTrialSub] + sub_tau[iSub]*(1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate threshold for high stakes
  mu_theta ~ normal(0, 1);          // theta: separate non-decision time for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);
	sd_theta ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);
	z_theta ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate threshold for high stakes
  real mu_theta;                                // theta: drift bonus for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate threshold for high stakes
  real<lower=0, upper=20> sd_theta;             // theta: drift bonus for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate threshold for high stakes
  real z_theta[nSub];                           // theta: drift bonus for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate threshold for high stakes
  real transf_mu_theta;                         // theta: drift bonus for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate threshold for high stakes
  real sub_theta[nSub];                         // theta: drift bonus for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  transf_mu_theta = mu_theta;
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: separate threshold for high stakes
    sub_theta[iSub] = mu_theta + z_theta[iSub] * sd_theta;                                  // theta: drift bonus for high stakes

    // Repeat for each trial of this subject:
    tau[firstTrial:lastTrial] = rep_vector(sub_tau[iSub],nTrial); 
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial alpha and tau:
      alpha[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_alpha[iSub]*(1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub])
      + sub_theta[iSub] * stakes[iTrialSub];


      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate threshold for high stakes
  mu_theta ~ normal(0, 1);          // theta: drift bonus for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);
	sd_theta ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);
	z_theta ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate non-decision time for high stakes
  real mu_theta;                                // theta: drift bonus for high stakes

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate non-decision time for high stakes
  real<lower=0, upper=20> sd_theta;             // theta: drift bonus for high stakes

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate non-decision time for high stakes
  real z_theta[nSub];                           // theta: drift bonus for high stakes
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate non-decision time for high stakes
  real transf_mu_theta;                         // theta: drift bonus for high stakes

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate non-decision time for high stakes
  real sub_theta[nSub];                         // theta: drift bonus for high stakes
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  transf_mu_theta = mu_theta;
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: separate non-decision time for high stakes
    sub_theta[iSub] = mu_theta + z_theta[iSub] * sd_theta;                                  // theta: drift bonus for high stakes

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub],nTrial); 
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub],nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial alpha and tau:
      tau[iTrialSub] = sub_pi[iSub]*stakes[iTrialSub] + sub_tau[iSub]*(1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub])
      + sub_theta[iSub] * stakes[iTrialSub];

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate non-decision time for high stakes
  mu_theta ~ normal(0, 1);          // theta: drift bonus for high stakes

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);
	sd_theta ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);
	z_theta ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
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
      return log1m((exp(-2 * alpha * beta * delta) - 1) / (exp(-2 * alpha * delta) - 1)); // analytical solution for integration (1-p because lower boundary)
    }
  }
}
data {

  int<lower=0> nStim;                     // total number of data points
  int<lower=0> nResp;                     // total number of data points
  int<lower=0> nTrial;                    // number of trials
  int<lower=0> nSub;                      // number of subject
  int<lower=0> nData;                     // total number of data points

  int<lower=1> stimuli[nData];            // cue identifier
  int<lower=0, upper=1> resp[nData];      // choices (1=upper bound, 0=lower bound)
  real<lower=0> rt[nData];                // RTs (in sec.)
  int<lower=-1, upper=1> outcome[nData];  // outcomes (-1 or 1)
  
  int<lower=0, upper=1> valence[nData];   // cue valence (1=Win, 0=Avoid)
  int<lower=0, upper=1> congruency[nData];// congruency (1=congruent, 0=incongruent)
  int<lower=0, upper=1> stakes[nData];    // stakes (1=high, 0=low)

  real RTbound;                           // lower bound or RT (e.g., 0.1 second; lower limit on tau)
  real minRT;                             // minimum RT of the observed data (upper limit for tau)
  matrix[nStim,nSub] Qi;                  // initial action values: for each cue (rows), subject (columns), used to initialize Q

}
parameters {    

  // Define hierarchical means (group-level).
  real mu_alpha;                                // alpha: boundary separation (threshold)
  real mu_tau;                                  // tau: nondecision time
  real mu_beta;                                 // beta: starting point bias
  real mu_deltaSlope;                           // deltaSlope: drift rate slope (scaling Q-value difference)
  real mu_deltaWin;                             // deltaWin: drift rate offset for Win cues
  real mu_deltaAvoid;                           // deltaAvoid: drift rate offset for Avoid cues
  real mu_epsilon;                              // epsilon: learning rate
  real mu_pi;                                   // pi: separate non-decision time for high stakes for congruent cues
  real mu_theta;                                // theta: separate non-decision time for high stakes for incongruent cues

  // Define hierarchical sds (group-level).
  real<lower=0, upper=20> sd_alpha;             // alpha: boundary separation (threshold)
  real<lower=0, upper=20> sd_tau;               // tau: nondecision time
  real<lower=0, upper=20> sd_beta;              // beta: starting point bias
  real<lower=0, upper=20> sd_deltaSlope;        // deltaSlope: drift rate slope (scaling Q-value difference)
  real<lower=0, upper=20> sd_deltaWin;          // deltaWin: drift rate offset for Win cues
  real<lower=0, upper=20> sd_deltaAvoid;        // deltaAvoid: drift rate offset for Avoid cues
  real<lower=0, upper=20> sd_epsilon;           // epsilon: learning rate
  real<lower=0, upper=20> sd_pi;                // pi: separate non-decision time for high stakes for congruent cues
  real<lower=0, upper=20> sd_theta;             // theta: separate non-decision time for high stakes for incongruent cues

  // Define subject-level parameters (standard-normal for non-centered parameterization):
  real z_alpha[nSub];                           // alpha: boundary separation (threshold)
  real z_tau[nSub];                             // tau: nondecision time
  real z_beta[nSub];                            // beta: starting point bias
  real z_deltaSlope[nSub];                      // deltaSlope: drift rate slope (scaling Q-value difference)
  real z_deltaWin[nSub];                        // deltaWin: drift rate offset for Win cues
  real z_deltaAvoid[nSub];                      // deltaAvoid: drift rate offset for Avoid cues
  real z_epsilon[nSub];                         // epsilon: learning rate
  real z_pi[nSub];                              // pi: separate non-decision time for high stakes congruent cues
  real z_theta[nSub];                           // theta: separate non-decision time for high stakes incongruent cues
  
} //end parameters
transformed parameters {

  // Define group-level means (transformed):
  real transf_mu_alpha;                         // alpha: boundary separation (threshold)
  real transf_mu_tau;                           // tau: nondecision time
  real transf_mu_beta;                          // beta: starting point bias
  real transf_mu_deltaSlope;                    // deltaSlope: drift rate slope (scaling Q-value difference)
  real transf_mu_deltaWin;                      // deltaWin: drift rate offset for Win cues
  real transf_mu_deltaAvoid;                    // deltaAvoid: drift rate offset for Avoid cues
  real transf_mu_epsilon;                       // epsilon: learning rate
  real transf_mu_pi;                            // pi: separate non-decision time for high stakes for congruent cues
  real transf_mu_theta;                         // theta: separate non-decision time for high stakes for incongruent cues

  // Define subject-level parameters (transformed):
  real <lower=0> sub_alpha[nSub];               // alpha: boundary separation (threshold)
  real <lower=0> sub_tau[nSub];                 // tau: nondecision time
  real <lower=0, upper=1> sub_beta[nSub];       // beta: starting point bias
  real <lower=0> sub_deltaSlope[nSub];          // deltaSlope: drift rate slope (scaling Q-value difference)
  real sub_deltaWin[nSub];                      // deltaWin: drift rate offset for Win cues
  real sub_deltaAvoid[nSub];                    // deltaAvoid: drift rate offset for Avoid cues
  real <lower=0, upper=1> sub_epsilon[nSub];    // epsilon: learning rate
  real <lower=0> sub_pi[nSub];                  // pi: separate non-decision time for high stakes for congruent cues
  real <lower=0> sub_theta[nSub];               // theta: separate non-decision time for high stakes for incongruent cues
  
  // Define trial-level parameters:
  vector[nData] alpha;                          // alpha: boundary separation (threshold)
  vector[nData] tau;                            // tau: nondecision time
  vector[nData] beta;                           // beta: starting point bias
  vector[nData] delta;                          // delta: drift rate (both positive and negative)

  // Define Q-values and PEs:
  matrix[nStim,nResp] Q;                        // Current action values for all cues: for each cue, for each response
  real PE;                                      // Prediction error on this trial

  // --------------------------------------------------
  
  // Transform group-level means for output:
  transf_mu_alpha = log1p_exp(mu_alpha);
  transf_mu_tau = log1p_exp(mu_tau);
  transf_mu_beta = inv_logit(mu_beta);
  transf_mu_deltaSlope = log1p_exp(mu_deltaSlope);
  transf_mu_deltaWin = mu_deltaWin;
  transf_mu_deltaAvoid = mu_deltaAvoid;
  transf_mu_epsilon = inv_logit(mu_epsilon);
  transf_mu_pi = log1p_exp(mu_pi);
  transf_mu_theta = log1p_exp(mu_theta);
  
  // --------------------------------------------------
  
  // Learning model:
  for (iSub in 1:nSub){                         // Loop over subjects
  
    int firstTrial = ((iSub-1)*nTrial+1); // first trial of this subject
    int lastTrial = (iSub*nTrial); // last trial of this subject

    // Transform parameters on subject-level:
    sub_alpha[iSub] = log1p_exp(mu_alpha + z_alpha[iSub] * sd_alpha);                       // alpha: transform to be positive
    sub_tau[iSub] = log1p_exp(mu_tau + z_tau[iSub] * sd_tau);                               // tau: transform to be positive
    sub_beta[iSub] = inv_logit(mu_beta + z_beta[iSub] * sd_beta);                           // beta: starting point bias bound to range 0-1
    sub_deltaSlope[iSub] = log1p_exp(mu_deltaSlope + z_deltaSlope[iSub] * sd_deltaSlope);   // deltaSlope: transform to be positive (don't invert QGo-QNoGo)
    sub_deltaWin[iSub] = mu_deltaWin + z_deltaWin[iSub] * sd_deltaWin;                      // deltaWin: drift rate offset for Win cues
    sub_deltaAvoid[iSub] = mu_deltaAvoid + z_deltaAvoid[iSub] * sd_deltaAvoid;              // deltaAvoid: drift rate offset for Avoid cues
    sub_epsilon[iSub] = inv_logit(mu_epsilon + z_epsilon[iSub] * sd_epsilon);               // epsilon: learning rate bound to range 0-1
    sub_pi[iSub] = log1p_exp(mu_pi + z_pi[iSub] * sd_pi);                                   // pi: separate non-decision time for high stakes for congruent cues
    sub_theta[iSub] = log1p_exp(mu_theta + z_theta[iSub] * sd_theta);                       // theta: separate non-decision time for high stakes for incongruent cues

    // Repeat for each trial of this subject:
    alpha[firstTrial:lastTrial] = rep_vector(sub_alpha[iSub], nTrial);
    beta[firstTrial:lastTrial] = rep_vector(sub_beta[iSub], nTrial); 

    // Initialize variables:
    Q[1:nStim,1] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    Q[1:nStim,2] = Qi[1:nStim,iSub];             // Initialize action values for all cues
    
    for(iTrial in 1:nTrial){                  // Loop over trials
      
      int iTrialSub = (iSub-1)*nTrial+iTrial; // index for this subject for this trial
      
      // Compute trial-by-trial tau given congruency and stakes:
      tau[iTrialSub] = (sub_pi[iSub] * congruency[iTrialSub] + sub_theta[iSub] * (1-congruency[iTrialSub])) * stakes[iTrialSub] + sub_tau[iSub] * (1-stakes[iTrialSub]);

      // Compute trial-by-trial drift rate:
      delta[iTrialSub] = sub_deltaSlope[iSub]*(Q[stimuli[iTrialSub],1] - Q[stimuli[iTrialSub],2]) 
      + sub_deltaWin[iSub] * valence[iTrialSub] 
      + sub_deltaAvoid[iSub] * (1 - valence[iTrialSub]);

      // Update Q-values:
      PE = outcome[iTrialSub] - Q[stimuli[iTrialSub],2-resp[iTrialSub]];                                        // Compute prediction error
      Q[stimuli[iTrialSub],2-resp[iTrialSub]] = Q[stimuli[iTrialSub],2-resp[iTrialSub]] + sub_epsilon[iSub]*PE; // Update Q-values
        
    } // end iTrial

  } // end iSub

} // end transformed parameters
model { // follow rather Kraemer et al.

  // Hierarchical parameters (group-level) priors of means:
  mu_alpha ~ normal(0, 1);          // alpha: boundary separation/threshold
  mu_tau ~ normal(0, 1);            // tau: nondecision time
  mu_beta ~ normal(0, 1);           // beta: starting point bias
  mu_deltaSlope ~ normal(5, 2);     // deltaSlope: drift rate slope (scaling Q-value difference)
  mu_deltaWin ~ normal(0, 1);       // deltaWin: drift rate offset for Win cues
  mu_deltaAvoid ~ normal(0, 1);     // deltaAvoid: drift rate offset for Avoid cues
  mu_epsilon ~ normal(0, 1);        // epsilon: learning rate
  mu_pi ~ normal(0, 1);             // pi: separate non-decision time for high stakes for congruent cues
  mu_theta ~ normal(0, 1);          // theta: separate non-decision time for high stakes for incongruent cues

  // Hierarchical parameters (group-level) priors of sds:
	sd_alpha ~ normal(0, 1);
	sd_tau ~ normal(0, 1);
	sd_beta ~ normal(0, 1);
	sd_deltaSlope ~ normal(0, 1);
	sd_deltaWin ~ normal(0, 1);
	sd_deltaAvoid ~ normal(0, 1);
	sd_epsilon ~ normal(0, 1);
	sd_pi ~ normal(0, 1);
	sd_theta ~ normal(0, 1);

  // Subject-level parameters in non-centered parameterization:
	z_alpha ~ normal(0, 1);
	z_tau ~ normal(0, 1);
	z_beta ~ normal(0, 1);
	z_deltaSlope ~ normal(0, 1);
	z_deltaWin ~ normal(0, 1);
	z_deltaAvoid ~ normal(0, 1);
	z_epsilon ~ normal(0, 1);
	z_pi ~ normal(0, 1);
	z_theta ~ normal(0, 1);

  // Loop over data, increment log-likelihood:
  for (iData in 1:nData) {
    target += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
  
} // end model
generated quantities {
  
  // Define log likelihood:
  real log_lik = 0;

  for (iData in 1:nData) {
    log_lik += wiener_gonogo_lpdf( rt[iData] | resp[iData], alpha[iData], tau[iData], beta[iData], delta[iData]);
  } // end iData
} // end generated quantities
"

# END
