#### MGNGUn Create Trial Sequence ####

rm(list = ls())
set.seed(19913010)
# setwd("//CNAS.RU.NL/U011144/Documents/AC_Teaching/OP3_SS_19/MGNG_Uncertainty/Task/TaskStimuli")
setwd("D:/Task_MGNGStakes/TrialHandler")
# Create practice trials
nTrials <- 5 # Very as function of trials in main task; one trial with invalid feedback added manually later 
nPart <- 4

# Initialize task conditions
allBlock <- 1:nPart
allCondition <- c(1,3,2,4)
allCue <- c("Pract_G2W","Pract_NG2W","Pract_G2A","Pract_NG2A")
allValence <- c(1,1,0,0)
allReqAction <- c(1,0,1,0)

for (iPart in 1:4){
  block <- rep(allBlock[iPart],nTrials)
  trialnr <- 1:nTrials
  stimulus <- rep(allCue[iPart],nTrials)
  condition <- rep(allCondition[iPart],nTrials)
  valence <- rep(allValence[iPart],nTrials)
  reqAction <- rep(allReqAction[iPart],nTrials)
  manipulation <- sample(c(rep(0,nTrials/2),rep(1,nTrials/2)),nTrials,replace = F)
  goValidity <- rep(1,nTrials)
  nogoValidity <- rep(1,nTrials)
  ISI <- rep(0.7,nTrials)
  ITI <- rep(1.5,nTrials)
  
  # Insert invalid cue
  iInvalid <- sample(2:nTrials,1)
  goValidity[iInvalid] <- 0
  nogoValidity[iInvalid] <- 0
  
  # Put together
  taskSettings = as.data.frame(cbind(block,trialnr,stimulus,condition,valence,reqAction,manipulation,goValidity,nogoValidity,ISI,ITI))
  write.csv(taskSettings,
            paste0("stimuluslist_pract_Part",iPart,".csv"), quote = F, row.names = F)
}

# # Read:
# setwd("D:/Task/TrialHandler")
# data <- read.csv("stimuluslist_test_stimNum_32_sub_001_Part1.csv")
# data

