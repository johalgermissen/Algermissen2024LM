# prepObs4stan.R
prepObs4stan <- function(inputDir, targetFile, filePattern = ".csv", sub2exclude = c()) {

  # inputDir     = rawdatadir
  # inputDir     = "/project/2420133.01/MGNGStakes/data/rawData/behavior/"
  # targetFile   = "/project/2420093.01/MGNGStakes/data/stan_results/data4stan.Rdata"
  # filePattern = ".csv"
  # sub2exclude = c(52)
  
  cat("Preprocess & save data to be used in STAN\n")
  
  # ============================================================================ #
  #### Load data: ####
  
  cat("Load data\n")
  data <- do.call("rbind", lapply(list.files(inputDir, pattern=filePattern, full=TRUE), read.csv, header=T))
  
  # ============================================================================ #
  #### Subject 40 missing, instead twice number 20 --> rename: ####
  
  data$Subject[data$Subject == 20 & data$Age == 23] <- 40
  
  ### Sort data:
  data <- data[order(data$Trialnr, decreasing = FALSE), ]
  data <- data[order(data$Subject, decreasing = FALSE), ]
  
  # ============================================================================ #
  #### Select variables: ####

  # names(data)
  data <- data[, c("Subject", "Trialnr", "Stimulus", "ReqAction", "Valence", "Manipulation", "Validity", "Response", "RT", "Outcome")]
  # drop Age, Gender, Hand, Condition
  # head(data)
  # table(data$Subject)
  
  # ============================================================================ #
  #### Turn stimulus identifier into numbers: ####
  
  cat("Turn stimulus identifier to numeric\n")
  data$Stimulus <- as.numeric(as.factor(data$Stimulus))
  # table(data$Stimulus)

  # ============================================================================ #
  #### Set lower RT bound, determine minimal RT for slowest subject: ####
  
  cat("Set lower RT bound, determine minimal RT for slowest subject\n")
  RTbound = 0 # minimum possible RT (assumption made)
  minRT = max(as.numeric(tapply(data$RT,data$Subject,min,na.rm = T))) # highest minimal RT of any participant observed
  
  # ============================================================================ #
  #### NoGo RTs to 0: ####
  
  cat("Set RTs for NoGo to 0 (stan does not accept NaN)\n")
  data$RT <- ifelse(is.na(data$RT), 0, data$RT)
  
  # library(lattice); densityplot(data$RT)
  stopifnot(sum(is.na(data$RT)) == 0)
  
  # ============================================================================ #
  #### Determine data dimension sizes: ####
  
  cat("Determine data dimensions\n")
  nStim <- length(unique(data$Stimulus)) # number cues
  nResp <- length(unique(data$Response)) # number responses
  nSub <- length(unique(data$Subject)) # number subjects
  nData <- nrow(data) # total number trials
  nTrial <- nData / nSub # number trials per subject
  nRep <- nTrial / nStim

  # ============================================================================ #
  #### Identify valence per cue, get initial Q values per subject: ####
  
  Qi <- matrix(NA, nStim, nSub)
  for (iSub in 1:nSub){ # iSub <- 1
    subRowIdx <- which(data$Subject == iSub) # identify rows of this subject
    Qi[1:nStim, iSub] <- as.numeric(tapply(data$Valence[subRowIdx], data$Stimulus[subRowIdx], mean)) - 0.5
  }

  # ============================================================================ #
  #### Reshape data into matrix format: ####
  
  cat("Extract from data frame\n")
  
  ## Extract into separate vectors:
  stimuli <- data$Stimulus # cue 
  reqAction <- data$ReqAction # required action
  valence <- data$Valence # valence
  congruency <- as.numeric(reqAction==valence)
  stakes <- data$Manipulation # stakes
  
  resp <- data$Response # response_all
  rt <- data$RT # RT
  
  validity <- data$Validity
  outcome <- data$Outcome # outcome

  # ============================================================================ #
  #### Remove subjects that should be excluded: ####
  
  cat("Remove excluded subjects from data \n")
  cat("Excluded subjects are:", sub2exclude, "\n")
  subIdx = rep(T, nSub) # repeat TRUE for each subject
  subIdx[sub2exclude] = F # set those to be excluded to FALSE
  rowIdx <- rep(subIdx, each = nTrial) # repeat TRUE/FALSE indices nTrial times
  
  # Update nSub, nData:
  nSub <- nSub - length(sub2exclude)
  nData <- nSub * nTrial
  
  subject <- rep(1:nSub, each = nTrial)
  trialnr <- rep(1:nTrial, times = nSub)
  
  ## Select rows for valid subjects:
  stimuli <- stimuli[rowIdx]
  reqAction <- reqAction[rowIdx]
  valence <- valence[rowIdx]
  congruency <- congruency[rowIdx]
  stakes <- stakes[rowIdx]
  
  resp <- resp[rowIdx]
  rt <- rt[rowIdx]
  
  validity <- validity[rowIdx]
  outcome <- outcome[rowIdx]

  Qi <- Qi[, subIdx]
  
  # ============================================================================ #
  #### Save as list:

  cat("Save data ...\n")
  dList <- list(nStim = nStim, nResp = nResp, nTrial = nTrial, nSub = nSub, nData = nData, 
                subject = subject, trialnr = trialnr,
                stimuli = stimuli, reqAction = reqAction, valence = valence, congruency = congruency, stakes = stakes,
                resp = resp, rt = rt, validity = validity, outcome = outcome, 
                Qi = Qi,
                RTbound = RTbound, minRT = minRT)

  save(dList, file = targetFile)
  
  cat("Saved data :-)\n")
} # end prepObs4stan.

# END OF FILE.