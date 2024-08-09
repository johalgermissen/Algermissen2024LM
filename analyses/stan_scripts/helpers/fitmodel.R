fitmodel <- function(model2fit, obs4stanFile, stanFitsFile, calcLogLik = TRUE, isTest = FALSE){
  #' Fit model in stan and save
  #' model2fit      - numeric, number of model to fit
  #' obs4stanFile   - string, complete path of data file to load
  #' stanFitsFile   - string, complete path of workspace to save all data objects to
  #' calcLogLik     - Boolean, output log_lik or not  
  #' No function output; will save relevant generated data objects to stanFitsFile
  #' Johannes Algermissen, 2023.
  
  require(stringr)
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  cat("Start fitting model",model2fit_str,"\n")
  
  ## Manually load data for testing:
  calcLogLik <- T
  
  # ======================================================================================================== #
  #### Load input data for Rstan: ####

  cat("Load input data\n")
  load(obs4stanFile) # contains most data
  
  # ======================================================================================================== #
  #### Specify which parameters STAN should give as output: #####
  
  cat("Select parameters for output\n")
  
  ## Different levels for output parameters:
  parLevel = c("mu_", "sd_", "transf_mu_", "sub_") # levels of output parameters
  
  ## All possible parameter types, sub-selection per model:
  if (suffix == "_standard"){
    parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "delta", "deltaWin", "deltaAvoid", "epsilon") 
    paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,9), c(1,2,4,5,6,9), c(1,2,3,6,7,8,9)) 
  } else if (suffix == "_deltaIntSlope"){
    parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta", "piCong", "piIncong") 
    paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                         c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                         c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12),
                         c(1,2,3,7,8,9,10,13,14)
                         )
  } else if (suffix == "_sepParam"){
    parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta") 
    paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                         c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                         c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), 
                         c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12)
                         )
  } else {
    stop("Output parameters for this suffix not defined")
  }
  cat(paste0("Parameters are ", paste0(parType[unlist(paramOut4Mod[model2fit])], collapse = ", "), "\n"))
  
  ## Put together all parameters:
  paramOutAll <- expand.grid(parType, parLevel) # first type, then level
  paramOutAll <- paramOutAll[, c(2, 1)] # invert columns
  paramOutAll <- do.call(paste0, paramOutAll) # concatenate
  paramOutAll <- c(paramOutAll, "log_lik") # add log likelihood
  
  ## Select parameters for selected model:
  nParType <- length(parType) # number of different parameter types (all possible parameters)
  paramIdx <- paramOut4Mod[[model2fit]] # extract indices for selected model
  paramOut <- paramOutAll[c(0*nParType+paramIdx, 1*nParType+paramIdx, 2*nParType+paramIdx, 3*nParType+paramIdx, 
                            length(paramOutAll))] # extract all parameter names for selected model
  
  paramInit = paramOut[1:(length(paramOut)-1)] # without log_lik (don't show log_lik for f1)
  if(!calcLogLik){paramOut = paramInit} # if calcLogLik set to FALSE in function call: don?t fit likelihoods (takes more time?)
  
  # ======================================================================================================== #
  #### Specify initial values: #####
  
  cat("Select initial values\n")
  init_fun <- eval(parse(text = paste0("init_fun_mod", model2fit_str)))
  
  # ======================================================================================================== #
  #### Print all relevant settings to console: ####
  
  # retrieve modelstring to use.
  mstr = get(paste("mstr", model2fit_str, sep = "")) # string to be used as input for stan fitting
  
  # some printing to check if thing go OK...
  cat("model2fit = ", model2fit_str, "\n")
  cat("Data file = ", obs4stanFile, "\n")
  cat("paramOut =" , paramOut, "\n")
  cat("STAN model string:\n")
  cat(mstr)
  
  # ======================================================================================================== #
  #### MODEL FITTING + PARAMETER ESTIMATION: ####
  
  cat(" ================================================================================================ \n")
  
  ## Specify seed for reproducibility of results.
  seed = 71 #  arbitrary, but setting a fixed seed makes models reproducible
  # in chains below in f2, seed is incremented by 1 for each subsequent chain, so the chains are actually different
  
  ## stan() requires the following input arguments:
  # model_code needs to be in string format for rstan
  # pars specifies which outputs you want returned: 
  # init = initial param values.  
  # If fit is not NA, the compiled model associated with the fitted result is re-used; thus the time that would otherwise be spent recompiling the C++ code for the model can be saved. So the 2step function here is to simply avoid recompiling the code 4 x when using 4 chains. 

  # ------------------------------------------------------------------------- #  
  ## 1) Run few iterations to compile the code for stan (checks already whether compiling works), so that later stages run much faster:
  # mstr = mstr14
  
  nIterTest = 10
  cat(paste0("Fit first initial model for compilation with only ", nIterTest, " iterations\n"))
  start.time <- Sys.time(); cat("Start at", as.character(start.time), "\n")
  f1 <- stan(model_code = mstr, data = dList, iter = nIterTest, chains = 1, init = init_fun,
             seed = seed, chain_id = 1, pars = paramInit)
  end.time <- Sys.time(); cat("Finish at", as.character(end.time), "\n")
  duration <- difftime(end.time, start.time); 
  cat(paste("\nFitting took", capture.output(duration)), "\n")
  options(max.print = 7200)
  print(f1)
  
  # ------------------------------------------------------------------------- #
  ## 2) Run actual model (4 long chains)

  if(isTest == FALSE){
    
    cat("Fit second, longer chains\n")

    ## Start taking time:
    start.time <- Sys.time()
  
    # ----------------------------------------------------------------------- #
    ## 1) Test chain with multiple cores:
    # f2 <- stan(fit = f1, data = dList, warmup = 200, iter = 600, chains = 4, cores = 4, init = init_fun,
    #            seed = seed+1, pars = paramOut) # fitting parameters for Model M8
    
    # ----------------------------------------------------------------------- #
    ## 2) Too short chain, will not converge for most models:
    # f2 <- stan(fit = f1, data = dList, warmup = 2000, iter = 6000, chains = 4, init = init_fun,
    #            seed = seed+1, pars = paramOut) # fitting parameters for Model M8
    #            # control = list(adapt_delta = 0.90)) # increase adapt_delta
    
    # ----------------------------------------------------------------------- #
    ## 3) Likely more reasonable length (only 2/3 of very long job) !!!!!!!!!!!!!!!
  
    # f2 <- stan(fit=f1, data=dList, warmup = 4000, iter = 8000, chains = 4, cores = 4, init = init_fun,
    #            seed=seed+1, pars=paramOut) # fitting parameters for all models by default
    f2 <- stan(fit=f1, data=dList, init = init_fun,
               # warmup = 500, iter = 1000, 
               # warmup = 4000, iter = 8000,
               warmup = 5000, iter = 10000,
               chains = 4, cores = 4,
               # control = list(adapt_delta = 0.90), # default 0.95
               # control = list(adapt_delta = .95, max_treedepth = 20),
               # control = list(adapt_delta = .95, max_treedepth = 15), # for M1-M8
               control = list(adapt_delta = .99, max_treedepth = 15), # for M6, 9-11
               seed = seed+1, pars = paramOut) # fitting parameters for all models by default
    # f2 <- stan(fit = f1, data = dList, warmup = 10000, iter = 20000, chains = 4, cores = 4, init = init_fun,
    #            seed = seed+1, pars = paramOut) # too long for study 2

    # ---------------------------------------------------------------------- #
    ## Stop taking time:
  
    end.time <- Sys.time(); cat("Finish at", as.character(end.time), "\n")
    duration <- difftime(end.time,start.time); cat(paste("\nFitting took", capture.output(duration)), "\n")
    options(max.print = 7200)
    print(f2)

  }
  
  # ------------------------------------------------------------------------- #
  ## 3) Save all relevant outputs as R workspace:

  if(isTest == FALSE){
    cat("Save data... \n")
    if(!exists("f2")){f2 <- f1} # if save, but f2 does not exist
    save(dList, mstr, f1, f2, duration, file = stanFitsFile)
    cat("Finished :-)\n")
  }
  
  } # end fitModel.
