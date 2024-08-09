# initVals4stan_sepParam.R
# Use separate parameters instead of difference parameters.
# Johannes Algermissen, 2023.

# Fontanesi et al. 2019 PBR: https://github.com/laurafontanesi/rlssm/blob/main/tests/notebooks/RLDDM_hierarchical_fitting.ipynb
# Kraemer et al. 2021 PBR: https://osf.io/4r25z/

nSub <- 54 # number parameters

cat(paste0("Initialize values for ", nSub, " subjects\n"))

alphaStart <- -0.18 # see Kraemer: used -0.18
tauStart <- -10 # Kraemer used -2.5 to achieve log(1 + exp(-2.5))) = 0.6 --> try more extreme
biasSDInit <- 0.1 # bias
driftSDInit <- 1 # drift intercepts
diffInit <- 0.1 # difference parameter
sdInit <- 0.001

# ======================================================================================================== #
#### init_fun_mod01: DDM (chance performance with bias) ####

init_fun_mod01 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha 	   = alphaStart, # used by Peter Kraemer, becomes around 0.60
    mu_tau   	   = tauStart, # used by Peter Kraemer, becomes around 0.08
    mu_beta  	   = rnorm(1, 0.5, biasSDInit), # use restricted sd
    mu_deltaInt  = rnorm(1, 0, driftSDInit),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha 	   = sdInit,
    sd_tau   	   = sdInit,
    sd_beta  	   = sdInit,
    sd_deltaInt  = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha  	   = rnorm(nSub, 0, 1),
    z_tau    	   = rnorm(nSub, 0, 1),
    z_beta   	   = rnorm(nSub, 0, 1),
    z_deltaInt   = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod02: RL-DDM #### 

init_fun_mod02 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaInt    = rnorm(1, 0, driftSDInit),
    mu_deltaSlope  = rnorm(1, 0, 1),
    mu_epsilon     = rnorm(1, 0, 1),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaInt    = sdInit,
    sd_deltaSlope  = sdInit,
    sd_epsilon     = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaInt     = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod03: RL-DDM with valence-dependent biases on starting point: #### 

init_fun_mod03 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_betaWin     = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_betaAvoid   = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaInt    = rnorm(1, 0, driftSDInit),
    mu_deltaSlope  = rnorm(1, 0, 1),
    mu_epsilon     = rnorm(1, 0, 1),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_betaWin     = sdInit,
    sd_betaAvoid   = sdInit,
    sd_deltaInt    = sdInit,
    sd_deltaSlope  = sdInit,
    sd_epsilon     = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_betaWin      = rnorm(nSub, 0, 1),
    z_betaAvoid    = rnorm(nSub, 0, 1),
    z_deltaInt     = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod04: RL-DDM with valence-dependent biases on drift: #### 

init_fun_mod04 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0, driftSDInit), # rnorm(1, 0, driftSDInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod05: RL-DDM with valence-dependent biases on drift and stakes on threshold: #### 

init_fun_mod05 <- function(chain_id = 1){
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0,regInit), # rnorm(1, 0,regInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = alphaStart,
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod06: RL-DDM with valence-dependent biases on drift and stakes on non-decision time: #### 

init_fun_mod06 <- function(chain_id = 1){
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0,regInit), # rnorm(1, 0,regInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = tauStart,
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod07: RL-DDM with valence-dependent biases on drift and stakes on bias: #### 

init_fun_mod07 <- function(chain_id = 1){
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0,regInit), # rnorm(1, 0,regInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod08: RL-DDM with valence-dependent biases on drift and stakes on drift: #### 

init_fun_mod08 <- function(chain_id = 1){
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0,regInit), # rnorm(1, 0,regInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = rnorm(1, 0, driftSDInit), 
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod09: RL-DDM with valence-dependent biases on drift and stakes on threshold & non-decision-time: #### 

## Stakes on threshold (pi) and non-decision time (theta)
init_fun_mod09 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0, driftSDInit), # rnorm(1, 0, driftSDInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = alphaStart,
    mu_theta       = tauStart,
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    sd_theta       = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1),
    z_theta        = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod10: RL-DDM with valence-dependent biases on drift and stakes on threshold & drift: #### 

## Stakes on threshold (pi) and drift (theta)
init_fun_mod10 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0, driftSDInit), # rnorm(1, 0, driftSDInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = alphaStart,
    mu_theta       = rnorm(1, 0, driftSDInit),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    sd_theta       = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1),
    z_theta        = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod11: RL-DDM with valence-dependent biases on drift and stakes on non-decision time & drift: #### 

## Stakes on non-decision time (pi) and drift (theta)
init_fun_mod11 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0, driftSDInit), # rnorm(1, 0, driftSDInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = tauStart,
    mu_theta       = rnorm(1, 0, driftSDInit),
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_pi          = sdInit,
    sd_theta       = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rnorm(nSub, 0, 1), # rep(0, nSub) # 
    z_theta        = rnorm(nSub, 0, 1)
  )
}

# ======================================================================================================== #
#### init_fun_mod12-14: RL-DDM with valence-dependent biases on drift and stakes effect on non-decision-time separately for congruent/incongruent cues: #### 

## Stakes on non-decision time with separate parameters for congruent and incongruent cues
init_fun_mod12 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(
    # Group-level means: fix tau to be small, alpha close to 1, beta close to 0.5
    mu_alpha       = alphaStart, # alphaStart used by Peter Kraemer, becomes around 0.60
    mu_tau         = tauStart, # -0.25 used by Peter Kraemer, becomes around 0.08 log(1+exp(mu_tau))
    mu_beta        = rnorm(1, 0.5, biasSDInit), # use restricted sd 
    mu_deltaSlope  = rnorm(1, 0, 1), # -3, # rnorm(1, 0, driftSDInit), # rnorm(1, 0, driftSDInit),
    mu_deltaWin    = rnorm(1, 0, driftSDInit),
    mu_deltaAvoid  = rnorm(1, 0, driftSDInit),
    mu_epsilon     = rnorm(1, 0, 1),
    mu_pi          = tauStart,
    mu_theta       = tauStart,
    # Group-level SDs: very small---subject level parameters largely ignored at the beginning
    sd_alpha       = sdInit,
    sd_tau         = sdInit,
    sd_beta        = sdInit,
    sd_deltaSlope  = sdInit,
    sd_deltaWin    = sdInit,
    sd_deltaAvoid  = sdInit,
    sd_epsilon     = sdInit,
    sd_piCong      = sdInit,
    sd_piIncong    = sdInit,
    # Subject level parameters: standard normally distributed:
    z_alpha        = rnorm(nSub, 0, 1),
    z_tau          = rnorm(nSub, 0, 1),
    z_beta         = rnorm(nSub, 0, 1),
    z_deltaSlope   = rnorm(nSub, 0, 1),
    z_deltaWin     = rnorm(nSub, 0, 1),
    z_deltaAvoid   = rnorm(nSub, 0, 1),
    z_epsilon      = rnorm(nSub, 0, 1),
    z_pi           = rep(0, nSub), #rnorm(nSub, 0, 1), #  
    z_theta        = rep(0, nSub) # rnorm(nSub, 0, 1) # 
  )
}

# END
