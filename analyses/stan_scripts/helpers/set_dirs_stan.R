## set_dirs.R
# Fish Task DISCOVERY scripts.
# Johannes Algermissen, 2023.


set_dirs_stan <- function(rootDir){
  
  ## Root directory:
  
  dirs <- list()
  
  # -------------------------------------------------------------------------- #
  ## Root directories:
  dirs$root       = "/project/2420133.01/MGNGStakes/"

  # -------------------------------------------------------------------------- #
  ## Scripts:
  dirs$scripts    = paste0(dirs$root, "2019_mgngstakes_johalg/stan_scripts/")
  dirs$models     = paste0(dirs$scripts, "models/")
  dirs$helpers    = paste0(dirs$scripts, "helpers/")
  dirs$simscripts = paste0(dirs$scripts, "simulations/")
  
  # -------------------------------------------------------------------------- #
  ## Data:
  dirs$data       = paste0(dirs$root, "data/")
  dirs$raw        = paste0(dirs$data, "rawData/behavior/")

  # -------------------------------------------------------------------------- #
  ## Results:
  dirs$results    = paste0(dirs$root, "results/stan_results/") 
  
  ## Plots:
  dirs$plots <- paste0(dirs$results, "plots/")
  dir.create(dirs$plots, showWarnings = FALSE)
  
  ## Simulations:
  dirs$sims       = paste0(dirs$results, "simulations/") 
  dir.create(dirs$sims, showWarnings = FALSE)

  ## Parameter recovery:
  dirs$paramRecov = paste0(dirs$results, "parameter_recovery/") 
  dir.create(dirs$paramRecov, showWarnings = FALSE)
  
  ## Parameter recovery:
  dirs$modelRecov = paste0(dirs$results, "model_recovery/") 
  dir.create(dirs$modelRecov, showWarnings = FALSE)
  
  return(dirs)
  
}
