#### MGNGUn Create Trial Sequence ####

rm(list = ls())
set.seed(19913010)
# setwd("//CNAS.RU.NL/U011144/Documents/AC_Teaching/OP3_SS_19/MGNG_Uncertainty/Task/TaskStimuli")
setwd("D:/Task_MGNGStakes/TrialHandler")

# -----------------------------------------------------------------------------------------------------
#### 1) Make cue sequence ####
makeCueSequence <- function(nPart,nStim,nRepStim,nTrialPart){
  RG <- array(NA, dim = c(nPart,nTrialPart))
  for (iPart in 1:nPart){
    # Randomize stimulus sequence: Continue shuffling until there are no more than 2 stimulus repetitions.
    stimseq <- rep(1:nStim,each=nRepStim) # repeat each stimulus as often as it should be repeated per Part
    stimRep <- 1 # just initialize stimRep to > 0
    while (length(stimRep)>0){
      tmpRG <- sample(stimseq) # shuffle vector
      # find trials where the stimulus is the same as on the next trial:
      # a) Use stimuli themselves (1-8):
#      repetitionIdx <- which(diff(tmpRG,lag=1)==0)
      # b) Compute cue conditions (1-4):
      condition <- (tmpRG-1) %% (nStim/2) + 1
      repetitionIdx <- which(diff(condition,lag=1)==0)
      # check if there are more than two stimulus repetitions:
      stimRep <- which(diff(repetitionIdx) == 1) # empty if no double repetitions
    } # end while loop
    RG[iPart,] <- tmpRG # save for each part
    rm(tmpRG,repetitionIdx,stimRep)
  }
  return(RG)
} # end MakeCueSequence
  
# -----------------------------------------------------------------------------------------------------
#### 2) Make Validity sequence ####
makeFBSequence <-function(nPart,nStim,nResp,nRepStim,prob){
  # Check how often valid and invalid feedback per Part:
  nPartValid = ceiling(nRepStim * prob) # number valid feedbacks for certain cue per Part
  if (!((nTrialPart * prob)%%1==0)){
    print(paste0('No integer number of valid responses; rounded to ',nPartValid/nRepStim*100,'% validity!'))
    print('No integer number of repetitions per cue per Part!')
  }
  nPartInvalid        <- nRepStim - nPartValid; # number invalid feedbacks per Part
  # ------------------------------------------------------------------------------
  fb <- array(NA, dim = c(nPart,nStim,nResp,nRepStim)) # Feedback
  for (iPart in 1:nPart){
    # 2) Randomize feedback sequence per cue and per response option: valid or invalid feedback for Go and NoGo
    feedseq <- c(rep(1,nPartValid),rep(0,nPartInvalid)) # initialize feedback sequence
    for (iStim in 1:nStim){ # sample for each stimulus
      for (iResp in 1:nResp){ # sample for each response
        fb[iPart,iStim,iResp,] <- sample(feedseq) # randomize sequence and store
      } # end iResp 
    } # iStim
  } # end iPart
  return(fb)
} # end makeFeedbackSequence

# -----------------------------------------------------------------------------------------------------
#### 3) Make ITI sequence ####
makeITISequence <-function(nPart,nStim,nRepStim, intercept, maxSteps){
  ITI <- array(NA, dim = c(nPart,nStim,nRepStim)) # ITI
  ITIvec <- c(0.1,-0.1,0.2,-0.2,0.3,-0.3,0.4,-0.4) 
  if (nRepStim %% 2 == 1){
    ITIsel <- c(0,ITIvec[1:min(nRepStim,maxSteps)]) # odd number of repetitions: add 0
  } else {
    ITIsel <- ITIvec[1:min(nRepStim,maxSteps)] # even number of repetitions
  }
  ITIvector <- rep(ITIsel,length.out=nRepStim)
  ITIintercept <- intercept
  for (iPart in 1:nPart){
    for (iStim in 1:nStim){ # sample for each stimulus
      # 3) Randomize ITI: gap between feedback and response.
      ITI[iPart,iStim,] <- ITIintercept + sample(ITIvector)
    } # end iStim
  } # end iPart
  return(ITI)
} # end makeFeedbackSequence
# -----------------------------------------------------------------------------------------------------
#### 4) Save files ####
saveSequence <- function(nPart,nStim,nTrialPart,nTrialTask,RG,fb,ISIAll,ITIAll,sub){
  block <- rep(1:nPart, each = nTrialPart)
  trialnr <- 1:nTrialTask # complete trial numbers
  # Put cue sequence together:
  cueAll = c(t(RG[1:nPart,])) # transpose and concatenate all parts
  condition <- (cueAll-1) %% (nStim/2) + 1
  # Retrieve valence, required action, manipulation:
  valence <-  ifelse(condition %in% c(1,3), 1, 0)
  reqAction <- ifelse(condition %in% c(1,2), 1, 0)
  manipulation <- ifelse(cueAll > nStim/2, 1, 0) # second half: manipulated 
  # Check that conditions are orthogonal:
#  rcor.test(cbind(valence,reqAction,manipulation)) # should be all 0 above diagonal
  if(sum(colSums(cor(cbind(valence,reqAction,manipulation)))) != 3){
    print("Task regressors not orthogonal!!")
  }
  # Initialize stimulus and validities:
  stimulus <- rep(NA,nTrialTask) # actual stimulus name, e.g. "A1"
  goValidity <- rep(NA,length(cueAll))
  nogoValidity <- rep(NA,length(cueAll))
  ITI <- rep(NA,length(cueAll))
  iPart <- 0
  # Retrieve stimulus, validity, iti for each trial:
  for (iTrial in 1:length(cueAll)){
    if(iTrial %in% c(seq(1,nTrialTask,nTrialPart))){ # for every start of new part:
      countStim <- rep(0,nStim) # initialize stimcount to zero
      iPart <- iPart + 1 # increment iPart
      stimMapping <- sample(unique(condition)) # random mapping of conditions on stimuli
    } # end iTrial
    stimSet <- LETTERS[iPart] # stimulus set to draw from (capital letters)
    stimulus[iTrial] <- paste0(stimSet,stimMapping[condition[iTrial]]) # concatenate set letter and stimulus number
    iCue <- cueAll[iTrial] # retrieve instance (complete information) of this stimulus
    countStim[iCue] <- countStim[iCue] + 1 # increase stimulus count
    goValidity[iTrial] <- fb[iPart,iCue,1,countStim[iCue]]
    nogoValidity[iTrial] <- fb[iPart,iCue,2,countStim[iCue]]
    ISI[iTrial] <- ISIAll[iPart,iCue,countStim[iCue]]
    ITI[iTrial] <- ITIAll[iPart,iCue,countStim[iCue]]
  } # end of iTrial
  # Concatenate all task settings into data frame:
  taskSettings = as.data.frame(cbind(block,trialnr,stimulus,condition,valence,reqAction,manipulation,goValidity,nogoValidity,ISI,ITI))
  
  # Save settings for each part separately:
  require(stringr)
  for (iPart in 1:nPart){
    write.csv(taskSettings[((iPart-1)*nTrialPart+1):(iPart*nTrialPart),],
              paste0("stimuluslist_test_stimNum_",nTrialTask,"_sub_",str_pad(sub, 3, pad = "0"),"_Part",iPart,".csv"), quote = F, row.names = F)
  } # end iPart
}
# -----------------------------------------------------------------------------------------------------
#### Task Settings ####

# Remember: cues with and without manipulation treated differently; so actually twice as many presentations!
nPart <- 4 # blocks, i.e. how often repeat task with separate cue sets
nStim <- 8 # number stimuli per task (i.e. cue conditions x manipulation conditions)
nResp <- 2 # number possible responses
prob <- .80 # feedback validity
# Use .80 for 5, 10, 15, 20
# Use 0.8571429 for 4, 8, 12

# Crucial setting to vary:
# Should be 10 for 20 repetitions per trial:
nRepStim <- 1 # number stimulus repetitions per Part (remember: effectively twice as much because separated for manipulation)
# Compute aggregate settings:
nTrialPart <- nStim * nRepStim # number stimuli per Part
nTrialTask <- nTrialPart*nPart  # number trials per entire task
nRepStimTotal <- nRepStim * nPart # total number of repetitions of certain cue per entire task

nSubs = 10
for (sub in 1:nSubs){
  print(paste0("Start subject ",sub))
  RG <- makeCueSequence(nPart,nStim,nRepStim,nTrialPart)
  fb <- makeFBSequence(nPart,nStim,nResp,nRepStim,prob)
  ISI <- makeITISequence(nPart,nStim,nRepStim,intercept=0.7,maxSteps=4)
  ITI <- makeITISequence(nPart,nStim,nRepStim,intercept=1.5,maxSteps=4)
  saveSequence(nPart,nStim,nTrialPart,nTrialTask,RG,fb,ISI,ITI,sub)
}

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# Read back in

setwd("D:/Task/TrialHandler")
data <- read.csv("stimuluslist_test_stimNum_256_sub_057_Part1.csv")
mean(data$ISI)
table(data$ISI)
mean(data$ITI)
table(data$ITI)
#data <- read.csv("stimuluslist_test_stimNum_320_sub_030_Part1.csv")
# head(data)
# -----------------------------------------------------------------------------------------------------
#### Plot sequence ####

RG <- makeCueSequence(nPart,nStim,nRepStim,nTrialPart)
fb <- makeFBSequence(nPart,nStim,nResp,nRepStim,prob)
ITI <- makeITISequence(nPart,nStim,nRepStim)

# Stimulus sequence per part:
for (iPart in 1:nPart){
  plot(RG[iPart,],type = "b",main = paste0("Stimulus sequence Part ",iPart))
  readline(prompt="Press [enter] to continue")
}

# Feedback validity sequency per stimulus per part:
for (iPart in 1:nPart){
  par(mar=c(1,1,1,1))
  par(mfrow=c(nStim,nResp))
  for (iStim in 1:nStim){
    for (iResp in 1:nResp){
      plot(fb[iPart,iStim,iResp,], type = "l", 
           main = paste0("Feedback sequence Part ", iPart, " Stimulus ", iStim, " Response", iResp))
    }
  }
  readline(prompt="Press [enter] to continue")
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow=c(1,1))
}

# fb[1,2,1,]

# END