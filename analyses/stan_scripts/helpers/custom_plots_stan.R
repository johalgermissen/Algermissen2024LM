#' custom_plots_stan.R
#' - Plot density of group-level parameters
#' - Plot violin plot of subject-level means

# ======================================================================================================== #
#### Retrieve names of all parameters for given model: ####

retrieve_parName <- function(model2fit, suffix){
  
  cat(paste0("Select all possible parameters for model ", model2fit, " suffix ", suffix, "\n"))
  
  if (suffix == "_standard"){
    if (model2fit == 1){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[")
    } else if(model2fit == 2){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "epsilon[")
    } else if(model2fit %in% c(3)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "delta[", "epsilon[")
    } else if(model2fit %in% c(4)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[")
    } else {
      stop("Unknown model2fit")
    }
  } else if (suffix == "_deltaIntSlope"){
    if (model2fit == 1){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[")
    } else if(model2fit == 2){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[", "deltaSlope[", "epsilon[")
    } else if(model2fit %in% c(3)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "deltaInt[", "deltaSlope[", "epsilon[")
    } else if(model2fit %in% c(4)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[")
    } else if(model2fit %in% c(5, 6, 7, 8)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[")
    } else if(model2fit %in% c(9, 10, 11)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[", "theta[")
    } else if(model2fit %in% c(12)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "piCong[", "piIncong[")
    } else {
      stop("Unknown model2fit")
    }
  } else if (suffix == "_sepParam"){
    if (model2fit == 1){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[")
    } else if(model2fit == 2){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "epsilon[")
    } else if(model2fit %in% c(3)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "delta[", "epsilon[")
    } else if(model2fit %in% c(4)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[")
    } else if(model2fit %in% c(5:8)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[")
    } else if(model2fit %in% c(9:14)){
      parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[", "theta[")
    } else {
      stop("Unknown model2fit")
    }
  } else {
    stop("Unknown suffix")
  }
  
  return(parNameVec)
  
}

# ======================================================================================================== #
#### Rename parameters: ####

rename_parameter <- function(inputString, model2fit, suffix = "_sepParam"){
  
  cat(paste0("Rename input string ", inputString, "\n"))
  if (suffix == "_sepParam"){
    if (model2fit == 5){
      inputString <- gsub("alpha", "alpha_Low", inputString)
      inputString <- gsub("pi", "alpha_High", inputString)
    } else if (model2fit == 6){
      inputString <- gsub("tau", "tau_Low", inputString)
      inputString <- gsub("pi", "tau_High", inputString)
    } else if (model2fit == 7){
      inputString <- gsub("beta", "beta_Low", inputString)
      inputString <- gsub("pi", "beta_High", inputString)
    } else if (model2fit == 8){
      inputString <- gsub("deltaSlope", "delta_Low", inputString)
      inputString <- gsub("pi", "delta_High", inputString)
    } else if (model2fit == 9){
      inputString <- gsub("alpha", "alpha_Low", inputString)
      inputString <- gsub("pi", "alpha_High", inputString)
      inputString <- gsub("tau", "tau_Low", inputString)
      inputString <- gsub("theta", "tau_High", inputString)
    } else if (model2fit == 10){
      inputString <- gsub("alpha", "alpha_Low", inputString)
      inputString <- gsub("pi", "alpha_High", inputString)
      inputString <- gsub("deltaSlope", "delta_Low", inputString)
      inputString <- gsub("theta", "delta_High", inputString)
    } else if (model2fit == 11){
      inputString <- gsub("tau", "tau_Low", inputString)
      inputString <- gsub("pi", "tau_High", inputString)
      inputString <- gsub("deltaSlope", "delta_Low", inputString)
      inputString <- gsub("theta", "delta_High", inputString)
    } else if (model2fit == 12){
      inputString <- gsub("tau", "tau_Low", inputString)
      inputString <- gsub("pi", "tau_HighCong", inputString)
      inputString <- gsub("theta", "tau_HighIncong", inputString)
    } else{
      cat("No renaming done\n")
    } 
  }
  
  cat(paste0("Output string is ", inputString, "\n"))
  return(inputString)
}

# ======================================================================================================== #
#### Add underscore after Greek letters if necessary, add backslashes, extract last 3-x components: ####

# inputString <- selParamVec[1]
# inputString <- parNameSel
# inputString <- parNameAxis
# string2Greek(inputString = selParamVec[1])

string2Greek <- function(inputString){
  
  cat(paste0("Input  string is ", inputString, "\n"))
  
  # --------------------------------------------------- #
  ## Add underscores after parameters:
  # cat("Add underscore after Greek letter if necessary\n")
  inputString <- gsub("deltaWin", "delta_Win", inputString)
  inputString <- gsub("deltaAvoid", "delta_Avoid", inputString)
  inputString <- gsub("deltaInt", "delta_Intercept", inputString)
  inputString <- gsub("deltaSlope", "delta_Slope", inputString)
  
  # --------------------------------------------------- #
  ## Add double backslash for Greek characters:
  # cat("Add 2 backslashes before Greek letters for rendering in LaTeX\n")
  inputString <- gsub("alpha", "\\\\alpha", inputString)
  inputString <- gsub("tau", "\\\\tau", inputString)
  inputString <- gsub("beta", "\\\\beta", inputString)
  inputString <- gsub("delta", "\\\\delta", inputString)
  inputString <- gsub("epsilon", "\\\\epsilon", inputString)
  inputString <- gsub("pi", "\\\\pi", inputString)
  inputString <- gsub("theta", "\\\\theta", inputString)
  
  # --------------------------------------------------- #
  ## Split string along underscores:
  # cat("Split string along underscores\n")
  inputString <- str_split(inputString, "_") # split along underscore
  inputString <- inputString[[1]] # extract from list
  
  ## Remove subject:
  nParam <- length(inputString)
  if (inputString[1] == "transf" & inputString[2] == "mu"){inputString <- inputString[3:nParam]}
  if (inputString[1] == "mu"){inputString <- inputString[2:nParam]}
  if (inputString[1] == "sd"){inputString <- inputString[2:nParam]}
  if (inputString[1] == "sub"){inputString <- inputString[2:nParam]}
  nParam <- length(inputString) # update
  
  ## Concatenate again:  
  if (nParam == 1){
    cat("Detected 1 component in inputString, add empty subscript\n")
    inputString <- paste0(inputString, "_{}")
  } else if (nParam == 2){
    cat("Detected 2 components in inputString\n")
    inputString <- paste(inputString[1], paste0("{", inputString[2], "}"), sep = "_")
    # } else  if (inputString[1] == "sub" & nParam == 3){
    #   cat("Detected 3 components, first one is sub, keep everything after sub\n")
    #   inputString <- paste(inputString[2], paste0("{", inputString[3], "}"), sep = "_")
    # } else  if (nParam == 2){ # if exactly 2: keep both
    #   cat("Detected 2 components, only keep last one\n")
    #   inputString <- inputString[nParam]
    # } else  if (nParam == 3){ # if exactly 3: keep only last (remove transf and mu)
    #     cat("Detected 3 components, only keep last one\n")
    #     inputString <- inputString[nParam]
  } else { # otherwise: remove transf and mu, keep 3rd, add 4th as subscript
    cat(paste0("InputString has ", nParam, " components: ", paste0(inputString, collapse = ", ", "\n")))
    stop("Detect > 2 components, unclear what to do")
    # cat("Detected > 3 components, treat 3rd as Greek, following ones as underscore\n")
    # subScript <- paste0(inputString[4:nParam], sep = "")
    # inputString <- paste(inputString[3], paste0("{", subScript, "}"), sep = "_")
  } # end if distinction
  
  cat(paste0("Output string is ", inputString, "\n"))
  return(inputString)
  
}
# ======================================================================================================== #
#### Determine y-axis limits given parameter name: ####

# inputString <- parNamePlot
# yLim <- find_yLim(inputString); yLim

find_yLim <- function(inputString){
  
  if (grepl("alpha", inputString)){
    # yLim <- c(1.10, 1.60)
    yLim <- c(1.20, 1.60)
  } else if (grepl("tau", inputString)){
    # yLim <- c(0.05, 0.40)
    yLim <- c(0.04, 0.40)
  } else if (grepl("beta", inputString)){
    # yLim <- c(0.15, 0.40)
    yLim <- c(0.14, 0.40)
  } else if (grepl("deltaSlope", inputString)){
    # yLim <- c(0.50, 4.50)
    yLim <- c(1.30, 13.30)
  } else if (grepl("delta", inputString)){
    # yLim <- c(0, 3.25)
    yLim <- c(0, 4.00)
  } else if (grepl("epsilon", inputString)){
    # yLim <- c(0.2, 0.8)
    yLim <- c(0.05, 0.60)
  } else {
    warning("Unknown yLim, return NA")
    yLim <- NA
  }
  
  cat(paste0("Assign yLim = c(", yLim[1], ", ", yLim[2], ")\n"))
  return(yLim)
  
}

# ======================================================================================================== #
#### Extract and plot group-level densities: ####

plot_group_density <- function(model2fit, suffix, pattern = "transf", isInteractive = TRUE, 
                               LWD = 4, FTS = 60, WIDTH = 640, HEIGHT = 640){
  # FTS formerly 70
  
  axisFTS <- 40
  
  require(ggplot2)
  require(latex2exp)
  
  require(plyr)
  require(stringr)

  # ----------------------------------------------------- #
  ## Create file name:
  
  cat("Create file name\n")
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  stanFitsFile <- paste0("M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
 
  # ----------------------------------------------------- #
  ## Load file:
  
  cat("Load file", stanFitsFile, " ...\n")
  load(paste0(dirs$results, stanFitsFile))
  posterior <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
  posterior <- as.data.frame(posterior)
  paramVec <- names(posterior)
  selParamVec <- paramVec[grep("transf", paramVec, fixed = T)] # sub-selection based on pattern
  
  # ----------------------------------------------------- #
  ## Select parameters of interest:
  
  selParamVec <- paramVec[grep(pattern, paramVec, fixed = T)] # sub-selection based on pattern
  nParam <- length(selParamVec)
  
  # ----------------------------------------------------- #
  ## Simple ggplot:

  for(iParam in 1:nParam){ # iParam <- 8 # for each parameter
    
    ## Name for plot:
    parNameSel <- selParamVec[iParam] # actual name of parameter
    posterior$x <- posterior[, parNameSel] # copy over
    cat("* ----------------------------------------------* \n")
    cat(paste0(">>> Start plotting ", parNameSel, "\n"))
    cat(paste0("xMin = ", round(min(posterior$x), 3), ", xMax = ", round(max(posterior$x), 3), "\n"))
    
    ## Name for saving plot:
    parNamePlot <- rename_parameter(parNameSel, model2fit, suffix) # rename
    
    ## Name for x-axis:
    parNameAxis <- rename_parameter(parNameSel, model2fit, suffix) # rename
    parNameAxis <- string2Greek(parNameAxis) # prepare for LaTeX printing
    
    ## Plot:
    p <- ggplot(posterior, aes(x = x)) +
      geom_density(color = "#C97106", fill = "#FF8C00", size = 4) +
      labs(x = TeX(paste0("$", parNameAxis, "$")), y = "", title = "") +
      theme_classic() + 
      theme(plot.margin = unit(c(0, 0.8, 0.5, 0.8), "cm")) +
      theme(axis.line = element_line(colour = 'black', size = 4),
            # axis.title.y = element_text(vjust = 1),
            axis.title = element_text(size = FTS), # +20 
            axis.text.x = element_text(colour = 'black', size = axisFTS, margin = margin(t = 10)),
            axis.text.y = element_text(colour = 'black', size = axisFTS, margin = margin(t = 10)),
            axis.ticks.length.x = unit(0.1, "cm"),
            # axis.text.x = element_text(vjust = -0.1),
            # axis.text.y = element_text(hjust = -0.1),
            title = element_text(size = FTS), # + 5
            legend.text = element_text(size = FTS))
    # print(p)
    
    ## Save:
    plotName <- paste0("density_ggplot_group_mod", model2fit_str, suffix, "_", parNamePlot, ".png")
    cat(paste0("Save as ", plotName, "\n"))
    png(paste0(dirs$plots, plotName),
        width = WIDTH, height = HEIGHT) #  type = "cairo")
    print(p)
    dev.off()
    
    ## Print again:
    print(p)
    if(isInteractive){
      readline(prompt = "Press [enter] to continue")
    }
  }
  cat("Finished all parameters :-)\n")
  # return(p)
}

# ======================================================================================================== #
#### Extract and plot subject-level means and densities: ####

plot_subject_means_density <- function(model2fit, suffix, isInteractive = TRUE,
                                       LWD = 4, FTS = 60, WIDTH = 640, HEIGHT = 640){

  ## FTS used to be 70
  
  require(ggplot2)
  require(latex2exp)
  
  require(plyr)
  require(stringr)
  
  # -------------------------------------------------------------------------- #
  ## Create file name:
  
  cat("Create file name\n")
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  stanFitsFile <- paste0("M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
  
  # -------------------------------------------------------------------------- #
  ## Load file:

  cat("Load file", stanFitsFile, " ...\n")
  load(paste0(dirs$results, stanFitsFile))
  
  # Process variables:
  posterior <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
  paramVec <- names(posterior) # extract names before creating data frame (not split per subject yet)
  selParamVec <- paramVec[grep("sub", paramVec, fixed = T)] # sub-selection based on pattern
  nParam <- length(selParamVec) # number parameters
  posterior <- as.data.frame(posterior) # convert to data frame
  
  # -------------------------------------------------------------------------- #
  ## Aggregate per subject:
  
  cat("Extract mean parameter per subject\n")
  nSub <- length(grep(selParamVec[1], names(posterior))) # number subjects
  
  ## Into data frame:
  subParam <- data.frame(subject = seq(1, nSub),
                         group = rep(1, nSub)) # initialize data frame
  
  ## Loop over parameters:
  for (iParam in 1:nParam){ # iParam <- 9
    
    ## Select parameter name:
    selPar <- selParamVec[iParam] # extract parameter name
    selParName <- str_split(selPar, "_") # decompose in substrings with _ as boundary
    selParName <- selParName[[1]] # unlist
    selParName <- paste(selParName[2:length(selParName)], collapse = "_") # extract reduced name
    
    selPar <- paste0(selPar, ".") # add dot to make sure that delta is not confused with deltaWin/deltaAvoid
    subNames <- names(posterior)[grep(selPar, names(posterior), fixed = T)] # # locate subject-specific variables
    subParam[, selParName] <- as.numeric(colMeans(posterior[, subNames])) # mean per column
  }
  
  # -------------------------------------------------------------------------- #
  ## Create plot using easyGgplot2:
  
  cat("Start plotting\n")
  subParamVec <- names(subParam)[3:ncol(subParam)] # extract variables to plot (skip subject and group)

  # LWD = 4; FTS = 70; WIDTH = 960; HEIGHT = 640;
  for(iParam in 1:nParam){ # iParam <- 5
    
    ## Name for plot:
    parNameSel <- subParamVec[iParam] # actual name of parameter
    subParam$y <- subParam[, parNameSel] # copy over
    # posterior$y <- posterior[, parNameSel] # copy over
    cat("* ----------------------------------------------* \n")
    cat(paste0(">>> Start plotting ", parNameSel, "\n"))
    cat(paste0("yMin = ", round(min(subParam$y), 3), ", yMax = ", round(max(subParam$y), 3), "\n"))
    
    ## Name for saving plot:
    parNamePlot <- rename_parameter(parNameSel, model2fit, suffix) # rename
    
    ## Name for x-axis:
    parNameAxis <- rename_parameter(parNameSel, model2fit, suffix) # rename
    parNameAxis <- string2Greek(parNameAxis) # prepare for LaTeX printing
    
    yLim <- find_yLim(parNamePlot)
    
    # ---------------------------------- #
    ## Plot:
    p <- ggplot(subParam, aes(x = group, y = y)) + 
      stat_ydensity(geom = "violin", fill = "#FEEFD9") + 
      geom_point(shape = 21, colour = "black", fill = "white", 
                 size = 5, stroke = 2, position=position_jitter(0.1)) +
      labs(x = TeX(paste0("$", parNameAxis, "$")),y = "Mean Est.",title = "") +
      coord_cartesian(ylim = yLim) + 
      theme_classic() + 
      theme(plot.margin = unit(c(0, 0.8, 0.5, 0.8), "cm")) +
      theme(axis.line = element_line(colour = 'black', size = LWD),
            # axis.title.y = element_text(vjust = 1),
            axis.title = element_text(size = FTS), #  + 20 
            axis.text = element_text(colour = 'black', size = FTS),
            axis.ticks.length.x = unit(0.2, "cm"),
            axis.text.x = element_blank(), # no x-axis text
            # axis.text.x = element_text(vjust = -0.1),
            # axis.text.y = element_text(hjust = -0.1),
            # axis.text.x = element_text(margin = margin(t = 10)),
            # axis.text.y = element_text(margin = margin(r = 10)),
            title = element_text(size = FTS), #  + 5
            legend.text = element_text(size = FTS))
    
    ## Save:
    png(paste0(dirs$plots, "density_ggplot_sub_mod", model2fit_str, suffix, "_", parNamePlot, ".png"), 
        width = WIDTH, height = HEIGHT) # type = "cairo")
    print(p)
    dev.off()
    
    ## Print again:
    print(p)
    if (isInteractive){
      readline(prompt = "Press [enter] to continue")
    }
  }
  cat("Finished all parameters :-)\n")
  # print(p)
  # return(p) 
}

# ======================================================================================================== #
#### Extract and plot all densities: ####

plot_stan_dens <- function(model2fit, suffix, WIDTH = 640, HEIGHT = 640){
  
  # -------------------------------------------------------------------------- #
  ## Create file name:
  
  cat("Create file name\n")
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  stanFitsFile <- paste0("M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
  
  # -------------------------------------------------------------------------- #
  ## Load file:
  
  cat("Load file", stanFitsFile, "\n")
  load(paste0(dirs$results, stanFitsFile))
  
  # -------------------------------------------------------------------------- #
  ## Select parameters:
  
  ## Extract parameter names for later:
  parAvgTable <- summary(f2)$summary
  paramVec <- rownames(parAvgTable) # extract names of all parameters
  
  parNameVec <- retrieve_parName(model2fit, suffix)
  
  cat(paste0("Parameters are ", paste(parNameVec, collapse = ", "), "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Loop and save:
  for (iPar in 1:length(parNameVec)){ # iPar <- 1
    parName <- parNameVec[iPar]
    parName <- sub("[", "", parName, fixed = TRUE)
    cat("* ----------------------------------------------* \n")
    cat(paste0("Start plotting ", parName, "\n"))
    selVec <- paramVec[grep(parName, paramVec, fixed = T)] # sub-selection based on pattern
    
    ## Create plot:
    p <- stan_dens(f2, pars = selVec)
    
    ## Save plot:
    png(paste0(dirs$plots, "density_mod", model2fit_str, suffix, "_", parName, ".png"),
        width = WIDTH, height = HEIGHT, type ="cairo")
    print(p)
    # stan_hist(f2, pars = "X[1]") # Histogram - rather use density
    dev.off()
    print(p)
  }
  cat("Finished all parameters :-)\n")
  # print(p)
  # return(p)
}

# ======================================================================================================== #
#### Extract and plot group-level densities: ####

plot_density_diff <- function(model2fit, suffix, var1, var2, LWD = 4, FTS = 60, WIDTH = 640, HEIGHT = 640){
  # FTS formerly 70
  
  require(ggplot2)
  require(latex2exp)
  
  require(plyr)
  require(stringr)
  
  # ----------------------------------------------------- #
  ## Create file name:
  
  cat("Create file name\n")
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  stanFitsFile <- paste0("M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
  
  # ----------------------------------------------------- #
  ## Load file:
  
  cat("Load file", stanFitsFile, " ...\n")
  load(paste0(dirs$results, stanFitsFile))
  posterior <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
  posterior <- as.data.frame(posterior)

  ## Create difference variable:
  posterior$diff <- posterior[, var1] - posterior[, var2]
  
  ## Make plot name:
  parNamePlot <- paste0(rename_parameter(var1, model2fit), "_-_", rename_parameter(var2, model2fit))
  
  ## Make axis name:
  parNameAxis <- paste0(string2Greek(rename_parameter(var1, model2fit)), " - ", string2Greek(rename_parameter(var2, model2fit)))

  # ----------------------------------------------------- #
  ## Simple ggplot:
  
  cat(paste0("Start plotting ", parNamePlot, " in raster plot (can be very slow) ...\n"))

  ## Plot:
  p <- ggplot(posterior, aes(x = diff)) +
    geom_density(color = "#C97106", fill = "#FF8C00", size = 4) +
    labs(x = TeX(paste0("$", parNameAxis, "$")), y = "", title = "") +
    theme_classic() + 
    theme(plot.margin = unit(c(0, 0.8, 0.5, 0.8), "cm")) +
    theme(axis.line = element_line(colour = 'black', size = 4),
          # axis.title.y = element_text(vjust = 1),
          axis.title = element_text(size = FTS), # +20 
          axis.text = element_text(colour = 'black', size = FTS),
          axis.ticks.length.x = unit(0.1, "cm"),
          # axis.text.x = element_text(vjust = -0.1),
          # axis.text.y = element_text(hjust = -0.1),
          axis.text.x = element_text(margin = margin(t = 10)),
          axis.text.y = element_text(margin = margin(r = 10)),
          title = element_text(size = FTS), # +5
          legend.text = element_text(size = FTS))
  
  ## Save:
  png(paste0(dirs$plots, "density_ggplot_diff_group_mod", model2fit_str, suffix, "_", parNamePlot, ".png"), 
      width = WIDTH, height = HEIGHT) # , type ="cairo")
  print(p)
  dev.off()
  
  ## Print again@
  print(p)
  cat("Finished plotting difference :-)\n")
  return(p)
  
} # end of function

# END OF FILE.