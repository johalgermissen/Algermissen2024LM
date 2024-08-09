## set_dirs.R
# Fish Task DISCOVERY scripts.
# Johannes Algermissen, 2023.


set_dirs <- function(rootDir){
  
  ## Root directory:
  
  dirs <- list()
  
  # -------------------------------------------------------------------------- #
  ## Root directories:
  
  dirs$root <- "/project/2420133.01/MGNGStakes/"

  # -------------------------------------------------------------------------- #
  ## Code:
  dirs$codeDir    <- paste0(dirs$root, "2019_mgngstakes_johalg/regression/")
  
  # -------------------------------------------------------------------------- #
  ## Data:
  
  ## Data:
  dirs$dataDir <- paste0(dirs$root, "data/")
  
  ## Raw data:
  dirs$rawDataDir <- paste0(dirs$dataDir, "rawData/behavior/")
  
  dirs$questDir <- paste0(dirs$dataDir, "rawData/questionnaires/")
  
  # -------------------------------------------------------------------------- #
  ## Data sets:
  dirs$processedDataDir <- paste0(dirs$dataDir, "processedData/")
  dir.create(dirs$processedDataDir, recursive = TRUE, showWarnings = FALSE) # recursive = TRUE)
  
  # -------------------------------------------------------------------------- #
  ## Results:
  dirs$resultDir <- paste0(dirs$root, "results/regression/")
  dir.create(dirs$resultDir, showWarnings = FALSE)
  
  ## Models:
  dirs$modelDir <- paste0(dirs$resultDir, "models/")
  dir.create(dirs$modelDir, showWarnings = FALSE)
  
  ## Plots:
  dirs$plotDir <- paste0(dirs$resultDir, "plots/")
  dir.create(dirs$plotDir, showWarnings = FALSE)
  
  return(dirs)
  
}