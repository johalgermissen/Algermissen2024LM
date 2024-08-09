#!/usr/bin/env Rscript
# ============================================================================ #
## 00_mgngstakes_functions_regression.R
## MGNGStakes functions for plotting and fitting regressions.
# Johannes Algermissen, 2023.

# ============================================================================ #
#### Read in behavioral data: #####

read_behavior <- function(){
  
  ## Find all raw data files:
  selPattern = ".csv"
  fileList <- list.files(dirs$rawDataDir, pattern = selPattern, full = TRUE)
  nFile <- length(fileList) # count subject
  if(nFile == 0){stop("No files found")}
  cat(paste0("Found ", nFile, " files\n"))
  
  ## Read data:
  data <- do.call("rbind", lapply(fileList, read.csv, header=T))
  
  ## Inspect:
  print(head(data, n = 10))
  print(table(data$Subject))
  cat("Finished! :-)\n")
  
  return(data)
}

# ============================================================================ #
#### Pre-process behavioral data: #####

wrapper_preprocessing <- function(data){
  #' Just a wrapper around various pre-processing functions specified below
  #' @param data    data frame with trial-level data
  #' @return data   same data frame with additional pre-processed variables
  
  ## Count variables:
  nVar <- ncol(data)
  
  # -------------------------------------------------------------------------- #
  ## Correct subject ID for double subject 20 (change one to subject ID 40):
  
  data$Subject[data$Subject == 20 & data$Age == 23] <- 40

  # -------------------------------------------------------------------------- #
  ## Sort data:
  
  data <- data[order(data$Trialnr, decreasing = FALSE), ]
  data <- data[order(data$Subject, decreasing = FALSE), ]
  
  # -------------------------------------------------------------------------- #
  ## Data markers:
  
  cat("Preprocess subject and session indices\n")
  
  data$subject_n <- data$Subject
  data$subject_f <- factor(data$subject_n)
  
  data$block_n <- data$Block
  data$block_f <- factor(data$block_n)
  
  data$trialnr_n <- data$Trialnr
  data$trialnr_f <- factor(data$trialnr_n)
  
  ## Extra grouping by block:
  data$subject_block_n <- data$subject_n*10 + data$block_n
  data$subject_block_f <- factor(data$subject_block_n)
  
  # -------------------------------------------------------------------------- #
  ### Stimuli:
  
  cat("Preprocess stimuli\n")
  
  ## Stimulus:
  data$stimulus_s <- data$Stimulus
  data$cue_n <- as.numeric(as.factor(data$Stimulus))
  data$cue_n <- data$cue_n - min(data$cue_n) + 1 # so minimum is 1
  data$cue_block_n <- (data$cue_n - 1) %% 4 + 1
  # table(data$cue_n, data$cue_block_n)
  
  ## Condition:
  data$condition_n <- data$Condition
  data$condition_f <- factor(data$condition_n, levels = 1:4, labels = c("G2W", "G2A", "NG2W", "NG2A"))
  data$condition_short1_f <- factor(data$condition_n, levels = 1:4, labels = c("G2W", "G2A", "N2W", "N2A"))
  data$condition_short2_f <- factor(data$condition_n, levels = 1:4, labels = c("G2W", "G2A", "N2W", "N2A"))
  
  ## Required action:
  data$reqAction_n <- data$ReqAction
  data$reqAction_f <- factor(data$reqAction_n, levels = c(1, 0), labels = c("Go", "NoGo"))

  ## Valence only:  
  data$valence_n <- data$Valence
  data$valence_f <- factor(data$valence_n, levels = c(1, 0), labels = c("Win", "Avoid"))
  data$valence_short_f <- factor(data$valence_n, levels = c(1, 0), labels = c("win", "avo"))
  
  ## Manipulation:
  data$stakes_n <- data$Manipulation
  data$stakes_f <- factor(data$stakes_n, levels = c(1, 0), labels = c("high", "low"))

  ## Combine conditions and stakes:
  data$condition_stakes_n <- data$condition_n*2 - data$stakes_n 
  data$condition_stakes_f <- factor(data$condition_stakes_n, levels = 1:8, 
                                    labels = c("G2Whigh", "G2Wlow", "G2Ahigh", "G2Alow", "NG2Whigh", "NG2Wlow", "NG2Ahigh", "NG2Alow"))
  data[1:50, c("condition_f", "stakes_f", "condition_stakes_f")]
  
  ## Validity:
  data$validity_n <- data$Validity
  data$validity_f <- factor(data$validity_n, levels = c(1, 0), labels = c("valid", "invalid"))
  
  # -------------------------------------------------------------------------- #
  ## Congruency valence and required response:
  data$reqCongruency_n <- ifelse(data$reqAction_n == data$valence_n, 1, 0)
  data$reqCongruency_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("congruent", "incongruent"))
  data$reqCongruency_short1_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("cong", "incong"))
  data$reqCongruency_short2_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("con", "inc"))
  
  ## Combine with stakes:
  data$reqCong_stakes_n <- (1 - data$reqCongruency_n) * 2 + (1 - data$stakes_n) + 1
  data$reqCong_stakes_f <- factor(data$reqCong_stakes_n, levels = c(1:4), labels = c("con_high", "con_low", "inc_high", "inc_low"))
  # data[1:20, c("reqCongruency_f", "stakes_f", "reqCong_stakes_f")]
  
  # -------------------------------------------------------------------------- #
  ### Responses:
  
  cat("Preprocess responses\n")
  
  data$response_n <- data$Response
  data$response_f <- factor(data$response_n, levels = c(1, 0), labels = c("Go", "NoGo"))

  # -------------------------------------------------------------------------- #
  ## Congruency valence and actual response:
  
  data$actCongruency_n <- ifelse(data$response_n == data$valence_n, 1, 0)
  data$actCongruency_f <- factor(data$actCongruency_n, levels = c(1, 0), labels = c("congruent", "incongruent"))
  data$actCongruency_short1_f <- factor(data$actCongruency_n, levels = c(1, 0), labels = c("con", "inc"))
  data$actCongruency_short2_f <- data$actCongruency_short1_f
  
  # -------------------------------------------------------------------------- #
  ## Accuracy:
  
  cat("Preprocess accuracy\n")
  
  data$ACC_n <- data$ACC
  data$ACC_f <- factor(data$ACC_n, levels = c(1, 0), labels = c("correct", "incorrect"))
  
  # -------------------------------------------------------------------------- #
  ## RTs:
  
  cat("Pre-process RTs\n")
  
  data$RT_n <- data$RT
  data$RT_log_n <- log(data$RT_n)
  
  data$RTcleaned_n <- data$RT
  data$RTcleaned_n[data$RT < 0.3] <- NA
  data$RTcleaned_log_n <- log(data$RTcleaned_n)
  
  # -------------------------------------------------------------------------- #
  ## RT in ms:
  
  data$RT_ms_n <- data$RT_n * 1000
  data$RT_ms_log_n <- log(data$RT_ms_n)
  
  # -------------------------------------------------------------------------- #
  ## Outcome:
  
  cat("Pre-process outcomes\n")
  
  data$outcome_n <- data$Outcome
  data$outcome_f <- factor(data$outcome_n, levels = c(1, 0, -1), labels = c("reward", "neutral", "punishment"))
  data$outcome_short_f <- factor(data$outcome_n, levels = c(1, 0, -1), labels = c("rew", "neu", "pun"))
  
  ## Relative outcomes (positive/ negative):
  data$outcome_rel_n <- ifelse(data$outcome_n == 1 | data$outcome_n == 0 & data$valence_n == 0, 1,
                               ifelse(data$outcome_n == -1 | data$outcome_n == 0 & data$valence_n == 1, 0, 
                                      NA))
  data$outcome_rel_f <- factor(data$outcome_rel_n, levels = c(1, 0), labels = c("positive", "negative"))
  data$outcome_rel_short_f <- factor(data$outcome_rel_n, levels = c(1, 0), labels = c("pos", "neg"))

  ## All outcomes (rewarded/ not rewarded/ not punished/ punished):
  data$outcome_all_n <- ifelse(data$outcome_n == 1, 1,
                               ifelse(data$outcome_n == 0 & data$valence_n == 1, 2,
                                      ifelse(data$outcome_n == 0 & data$valence_n == 0, 3,
                                             ifelse(data$outcome_n == -1, 4,
                                                    NA))))
  data$outcome_all_f <- factor(data$outcome_rel_n, levels = 1:4, labels = c("rewarded", "not rewarded", "not punished", "punished"))
  data$outcome_all_short_f <- factor(data$outcome_rel_n, levels = 1:4, labels = c("rew", "¬rew", "¬pun", "pun"))
  
  # -------------------------------------------------------------------------- #
  ## Trial-number within block, cue repetition:
  
  data <- add_cueRep(data)
  data$cueRep_f <- factor(data$cueRep_n)
  
  ## Add factors:
  data$response_last_f <- factor(data$response_last_n, levels = c(1, 0), labels = c("Go", "NoGo"))
  data$outcome_last_f <- factor(data$outcome_last_n, levels = c(1, 0, -1), labels = c("reward", "neutral", "punishment"))
  data$outcome_last_short_f <- factor(data$outcome_last_n, levels = c(1, 0, -1), labels = c("rew", "neu", "pun"))
  # View(data[, c("Subject", "Block", "Trialnr", "trialnr_block_n", "Stimulus", "cue_block_n", "cueRep_n")])
  
  ## Relative outcomes last trial (positive/ negative):
  data$outcome_last_rel_n <- ifelse(data$outcome_last_n == 1 | data$outcome_last_n == 0 & data$valence_n == 0, 1,
                               ifelse(data$outcome_last_n == -1 | data$outcome_last_n == 0 & data$valence_n == 1, 0, 
                                      NA))
  data$outcome_last_rel_f <- factor(data$outcome_last_rel_n, levels = c(1, 0), labels = c("positive", "negative"))
  data$outcome_last_rel_short_f <- factor(data$outcome_last_rel_n, levels = c(1, 0), labels = c("pos", "neg"))
  
  # -------------------------------------------------------------------------- #
  ## Response last trial:
  
  cat("Create response last trial\n")
  
  data <- add_lag(data, "response_n")
  data$response_lag1_f <- factor(data$response_lag1_n, levels = c(1, 0), labels = c("Go", "NoGo"))
  # data[, c("trialnr_block_n", "response_n", "response_lag1_n")]

  # -------------------------------------------------------------------------- #
  ## Repeat/switch:
  
  cat("Create repeat/ switch variable\n")
  
  data$repeat_n <- ifelse(data$response_n == data$response_last_n, 1, 
                          ifelse(data$response_n != data$response_last_n, 0, NA))
  data$repeat_f <- factor(data$repeat_n, levels = c(1, 0), labels = c("repeat", "switch"))
  
  # -------------------------------------------------------------------------- #
  ## Outcome last trial:
  
  cat("Create outcome last trial\n")
  
  data <- add_lag(data, "outcome_n")
  data$outcome_lag1_f <- factor(data$outcome_lag1_n, levels = c(1, 0), labels = c("reward", "punishment"))
  data$outcome_lag1_short_f <- factor(data$outcome_lag1_n, levels = c(1, 0), labels = c("rew", "pun"))
  # data[, c("trialnr_block_n", "outcome_n", "outcome_lag1_n")]
  
  # -------------------------------------------------------------------------- #
  ## Remove old variables:
  data <- data[, (nVar + 1):ncol(data)]
  
  cat("Finished! :-)\n")
  return(data)
}

# ============================================================================ #
#### Select data, standardize variables: ####

select_standardize <- function(data, sub2excl = c(52)){
  #' Exclude subjects, standardize numerical variables of remaining data.
  #' @param data      data frame with trial-level data
  #' @param sub2excl  vector of integers, IDs of subjects to exclude (default: none).
  #' @return modData  data frame selected subjects removed and variables standardized. 
  
  # -------------------------------------------------------------------------- #
  ### Select data:
  
  if (length(sub2excl) == 0){
    cat("Retain data of all subjects\n")
  } else {
    cat(paste0("Select data, exclude subjects ", paste0(sub2excl, collapse = ", "), "\n"))
  }
  ## Select subject:
  modData <- subset(data, !(subject_n %in% sub2excl))
  
  # -------------------------------------------------------------------------- #
  ## Standardize variables:
  
  cat("Standardize relevant numeric variables...\n")
  
  modData$trialnr_z <- as.numeric(scale(modData$trialnr_n))
  modData$trialnr_block_z <- as.numeric(scale(modData$trialnr_block_n))
  modData$cueRep_z <- as.numeric(scale(modData$cueRep_n))
  
  modData$RT_z <- as.numeric(scale(modData$RT_n))
  modData$RT_log_z <- as.numeric(scale(modData$RT_log_n))
  modData$RT_ms_log_z <- as.numeric(scale(modData$RT_ms_log_n))
  modData$RTcleaned_z <- as.numeric(scale(modData$RTcleaned_n))
  modData$RTcleaned_log_z <- as.numeric(scale(modData$RTcleaned_log_n))
  
  cat("Finished! :-)\n")
  
  return(modData)
  
}

# ============================================================================ #
#### Load and pre-process questionnaire data: #####

load_preprocess_questionnaires <- function(){
  
  ## Load data:
  cat("Load questionnaire data\n")
  fileName <- "results-survey367265.csv"
  data <- read.csv(paste0(dirs$questDir, fileName))

  ## Information on questions asked:  
  questionVec <- names(data)
  
  ## Assign names:
  cat("Rename raw data\n")
  names(data) <- c("RespID", "Date", "LastPage", "StartLanguage", "subject_n", "age_f", "gender_g", "forward_span_n", "backward_span_n",
                   "BIS_01", "BIS_02", "BIS_03", "BIS_04", "BIS_05", "BIS_06", "BIS_07", "BIS_08", "BIS_09", "BIS_10", "BIS_11",
                   "neuroticism_01", "neuroticism_02", "neuroticism_03", "neuroticism_04", "neuroticism_05", "neuroticism_06", "neuroticism_07", "neuroticism_08", "neuroticism_09", "neuroticism_10",
                   "neuroticism_11", "neuroticism_12", "neuroticism_13", "neuroticism_14", "neuroticism_15", "neuroticism_16", "neuroticism_17", "neuroticism_18", "neuroticism_19", "neuroticism_20",
                   "hypothesis_s", "strategies_indicated_n", "strategies_text_s", "difficulty_indicated_n", "difficulty_which_s", "difficulty_text_s",
                   "Comments", "Debriefing")
  
  ### Filter out NAs and too high subject numbers (tests)
  validSubs <- which(!(is.na(data$subject_n)) & data$subject_n %in% 1:55)
  data <- data[validSubs,]
  
  ### Check subject numbers:
  # length(data$subject_n) # 56 subjects -- one too much
  # length(unique(data$subject_n)) # only 54 unique subjects, so 2 subject numbers twice?
  # which(table(data$subject_n) != 1) # 25 and 47 double
  
  ### Correct subject numbers:
  ## Subject 25 twice: Second entry of subject 25 (row 26) is subject 26
  cat("Relable second instance of subject 25 as subject 26\n")
  doubleSubs <- which(data$subject_n == 25)
  data$subject_n[doubleSubs[2]] <- doubleSubs[2]
  
  ## Subject 47 twice:
  # data[which(data$subject_n == 47), ]
  # 47 1st entry only filled until Difficulty_Explanation, missed comments
  # --> Keep 1st entry, manually add comments from 2nd entry, delete 2nd entry
  cat("Combine entries for subject 47\n")
  doubleSubs <- which(data$subject_n == 47)
  data$Comments[doubleSubs[1]] <- data$Comments[doubleSubs[2]]
  data <- data[-doubleSubs[2], ]
  
  ## Final check:
  stopifnot(length(data$subject_n) == 55)
  stopifnot(length(unique(data$subject_n)) == 55)
  stopifnot(length(which(table(data$subject_n) != 1)) == 0)
  # length(data$subject_n) # 55
  # length(unique(data$subject_n)) # 55
  # which(table(data$subject_n) != 1) # none double
  
  # -------------------------------------------------------------------------- #
  ### Compute mean span:
  
  cat("Compute mean memory span\n")
  data$mean_span_n <- (data$forward_span_n + data$backward_span_n)/2
  
  # -------------------------------------------------------------------------- #
  ### Compute mean of Barratt Impulsiveness Scale ("non-planning" subscale, 11 items):
  
  # questionVec[c(10:20)]
  # items 7, 9, 10 framed as default
  # items 1, 2, 3, 4, 5, 6, 8, 11 reversed
  cat("Compute mean of Barratt Impulsiveness scale (non-planning subscale)\n")
  data$BIS_total_n <- (6 - data$BIS_01 + 6 - data$BIS_02 + 6 - data$BIS_03 + 6 - data$BIS_04 + 6 - data$BIS_05 +  
                       6 - data$BIS_06 + data$BIS_07 + 6 - data$BIS_08 + data$BIS_09 + data$BIS_10 + 6 - data$BIS_11)/11 
  
  # -------------------------------------------------------------------------- #
  ### Compute mean of Neuroticism subscale (20 items):
  
  # questionVec[c(21:40)]
  # items 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16 framed as default
  # items 7, 8, 9, 10, 17, 18, 19, 20 reversed
  cat("Compute mean of neuroticism scale\n")
  data$neuroticism_total_n <- (data$neuroticism_01 + data$neuroticism_02 + data$neuroticism_03 + data$neuroticism_04 + data$neuroticism_05 +
                                 data$neuroticism_06 + 6 - data$neuroticism_07 + 6 - data$neuroticism_08 + 6 - data$neuroticism_09 + 6 - data$neuroticism_10 +
                                 data$neuroticism_11 + data$neuroticism_12 + data$neuroticism_13 + data$neuroticism_14 + data$neuroticism_15 +
                                 data$neuroticism_16 + 6 - data$neuroticism_17 + 6 - data$neuroticism_18 + 6 - data$neuroticism_19 + 6 - data$neuroticism_20)/20
  
  return(data)
  }

  
# ============================================================================ #
#### Retrieve default values for plots: #####

retrieve_plot_defaults <- function(input){
  #' Retrieve default for given plotting parameter.
  #' @param input scalar string, name of parameter for which to set default.
  #' @return scalar numeric, default value of that parameter
  # retrieve_plot_defaults("FTS")
  
  if (input == "FTS"){
    output <- 28
  } else if (input == "LWD"){
    output <- 1.5
  } else if (input == "dotSize"){
    output <- 0.5
  } else if (input == "dodgeVal"){
    output <- 0.6
  } else {
    stop("Unknown input to function retrieve_plot_defaults()")
  }
  
  return(output)
}

# ============================================================================ #
#### Automatically convert variable names to pretty names for plots: #####

substitute_label <- function(labels){
  #' Substitute certain manually defined factor names for prettier names.
  #' @param labels vector of strings with names of factors in model.
  #' @return same vector with certain strings substituted.
  
  for (iItem in 1:length(labels)){
    
    cat(paste0("Automatically substitute factor labels for ", labels[iItem], " according to manual mapping\n"))
    
    labels[iItem] <- gsub("(Intercept)", "Intercept", labels[iItem])
    labels[iItem] <- gsub("X.Intercept.", "Intercept", labels[iItem])
    
    #### Indices:
    #### Indices:
    labels[iItem] <- gsub("subject_n", "Subject", labels[iItem])
    labels[iItem] <- gsub("subject_f", "Subject", labels[iItem])
    
    labels[iItem] <- gsub("trialnr_block_n", "Trial number within block", labels[iItem])
    labels[iItem] <- gsub("trialnr_block_f", "Trial number within block", labels[iItem])
    
    labels[iItem] <- gsub("trialnr_n", "Trial number", labels[iItem])
    labels[iItem] <- gsub("trialnr_f", "Trial number", labels[iItem])
    
    labels[iItem] <- gsub("block_n", "Block number", labels[iItem])
    labels[iItem] <- gsub("block_f", "Block number", labels[iItem])
    
    labels[iItem] <- gsub("cueRep_n", "Cue repetition", labels[iItem])
    labels[iItem] <- gsub("cueRep_f", "Cue repetition", labels[iItem])    
    ### Stimuli:
    labels[iItem] <- gsub("valence_n", "Valence", labels[iItem])
    labels[iItem] <- gsub("valence_f", "Valence", labels[iItem])
    labels[iItem] <- gsub("valence_short_f", "Val.", labels[iItem])
    
    labels[iItem] <- gsub("reqAction_n", "Req. action", labels[iItem])
    labels[iItem] <- gsub("reqAction_f", "Req. action", labels[iItem])
    labels[iItem] <- gsub("reqAction_short_f", "Req. act.", labels[iItem])

    labels[iItem] <- gsub("reqCong_stakes_n", "Congruency x Stakes", labels[iItem])
    labels[iItem] <- gsub("reqCong_stakes_f", "Congruency x Stakes", labels[iItem])
    
    labels[iItem] <- gsub("stakes_n", "Stakes", labels[iItem])
    labels[iItem] <- gsub("stakes_f", "Stakes", labels[iItem])
    
    labels[iItem] <- gsub("reqCongruency_n", "Congruency", labels[iItem])
    labels[iItem] <- gsub("reqCongruency_f", "Congruency", labels[iItem])
    labels[iItem] <- gsub("reqCongruency_short1_f", "Congruency", labels[iItem])
    labels[iItem] <- gsub("reqCongruency_short2_f", "Cong.", labels[iItem])
    
    labels[iItem] <- gsub("actCongruency_n", "Congruency", labels[iItem])
    labels[iItem] <- gsub("actCongruency_f", "Congruency", labels[iItem])
    labels[iItem] <- gsub("actCongruency_short_f", "Cong.", labels[iItem])

    labels[iItem] <- gsub("condition_n", "Cue condition", labels[iItem])
    labels[iItem] <- gsub("condition_f", "Cue condition", labels[iItem])
    labels[iItem] <- gsub("condition_short1_f", "Cue condition", labels[iItem])
    labels[iItem] <- gsub("condition_short2_f", "Cue \ncond.", labels[iItem])
    
    ### Responses:
    labels[iItem] <- gsub("simResp_n", "Sim. p(Go)", labels[iItem])
    labels[iItem] <- gsub("simResp_cleaned_n", "Sim. p(Go)", labels[iItem])
    labels[iItem] <- gsub("response_n", "p(Go)", labels[iItem])
    labels[iItem] <- gsub("response_f", "Response", labels[iItem])

    labels[iItem] <- gsub("simACC_n", "sim. p(Correct)", labels[iItem])
    labels[iItem] <- gsub("ACC_n", "p(Correct)", labels[iItem])
    labels[iItem] <- gsub("ACC_f", "Accuracy", labels[iItem])
    
    labels[iItem] <- gsub("repeat_n", "p(repeat last response)", labels[iItem])
    labels[iItem] <- gsub("repeat_f", "Last. response", labels[iItem])
    labels[iItem] <- gsub("repeat_short_f", "Last. resp.", labels[iItem])

    labels[iItem] <- gsub("simRT_n", "Sim. RT (in sec.)", labels[iItem])
    labels[iItem] <- gsub("simRT_cleaned_n", "Sim. RT (in sec.)", labels[iItem])
    
    labels[iItem] <- gsub("RT_n", "RT (in sec.)", labels[iItem])
    labels[iItem] <- gsub("RTcleaned_n", "RT (in sec.)", labels[iItem])
    labels[iItem] <- gsub("RTcleaned_z", "RT (z)", labels[iItem])
    
    labels[iItem] <- gsub("RT_ms_n", "RT (in ms)", labels[iItem])
    labels[iItem] <- gsub("RT_log_n", "RT (log(sec.))", labels[iItem])
    labels[iItem] <- gsub("RT_ms_log_n", "RT (log(ms))", labels[iItem])
    
    ### Outcome:
    labels[iItem] <- gsub("outcome_n", "Outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_f", "Outcome", labels[iItem])
    
    labels[iItem] <- gsub("outcome_lag1_n", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_lag1_f", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_lag1_short_f", "Prev. out.", labels[iItem])
    
    labels[iItem] <- gsub("outcome_last_n", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_last_f", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_last_short_f", "Prev. out.", labels[iItem])
    labels[iItem] <- gsub("outcome_last_rel_n", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_last_rel_f", "Previous outcome", labels[iItem])
    labels[iItem] <- gsub("outcome_last_rel_short_f", "Prev. out.", labels[iItem])
    
    ## General replacements:
    labels[iItem] <- gsub("1", "", labels[iItem])
    labels[iItem] <- gsub(":", "\nx ", labels[iItem])
    
  }
  
  return(labels)
  
}

# ============================================================================ #
#### Retrieve colour given independent variable: #####

retrieve_colour <- function(input){
  #' Retrieve colour scheme given variable name (pattern) as input.
  #' If none found, then use colorbrewer to create unidimensional YlOrRd color scheme,
  #' which is repeated as necessary to achieve number of required levels.
  #' @param input scalar string, pattern within variable name to match to colour scheme.
  #' @return output vector of strings, colours to use.
  
  if (grepl("reqAction", input, fixed = TRUE)){ 
    output <- c("#B2182B", "#2166AC")
  } else if (grepl("valence", input, fixed = TRUE)){ 
    # output <- c("#76F013", "#ED462F") # old rich/poor colors
    output <- c("#007174", "#FF654E") # checked for red/green-friendliness
  } else if (grepl("condition", input, fixed = TRUE)){ 
    # output <- c("#76F013", "#ED462F", "#76F013", "#ED462F") # old rich/poor colors
    output <- c("#007174", "#FF654E", "#007174", "#FF654E") # checked for red/green-friendliness
  } else if (grepl("reqCongruency", input, fixed = TRUE)){ 
    output <- c("#4DAC26", "#D01C8B")
  } else if (grepl("actCongruency", input, fixed = TRUE)){ 
    output <- c("#4D9221", "#C51B7D")
  } else if (grepl("reqCong_stakes", input, fixed = TRUE)){ 
    output <- c("#574571", "#DEC5DA", "#574571", "#DEC5DA")
  } else if (grepl("stakes", input, fixed = TRUE)){ 
    # output <- c("orange", "purple")
    output <- c("#574571", "#DEC5DA") # dark purple, light puprle

  } else  if (grepl("response", input, fixed = TRUE)){ 
    output <- c("#B2182B", "#2166AC")
  } else if (grepl("ACC", input, fixed = TRUE)){ 
    output <- c("#721B3E", "#FFB2F3") # dark red/ light red from MetBrewer: Austria(7)
  } else  if (grepl("repeat", input, fixed = TRUE)){ 
    output <- c("#453D8D", "#B0B0FF") # dark violett/ light violett from MetBrewer: Redon(12) #59385C #A1A1FF #AB84A5
  } else  if (grepl("RTcleaned", input, fixed = TRUE)){ 
    output <- c("#183571", "#7EC5F4") # dark blue/ light blue from MetBrewer: Manet(11) #7EC5F4 #8CC8BC
    
  } else if (grepl("outcome_rel", input, fixed = TRUE)){ 
    output <- c("#009933", "#CC0000")
  } else if (grepl("outcome_last_rel", input, fixed = TRUE)){ 
    output <- c("#009933", "#CC0000")
  } else if (grepl("outcome", input, fixed = TRUE)){ 
    output <- c("#009933", "grey", "#CC0000")
  } else {
    if (exists("plotData")){ # check if plotData exists
      
      require(RColorBrewer)
      cat("Could not find variable name, but found plotData, retrieve number of variable levels\n")
      
      ## Count labels of variable:
      nLevel <- length(unique(plotData[, input]))
      cat(paste0("Found ", nLevel, " levels of variable ", input, " in plotData\n"))
      
      ## Retrieve number colours and repetitions:      
      cmap <- "YlOrRd"
      nColour <- min(nLevel, 9)
      nRep <- ceil(nLevel/9)
      cat(paste0("Retrieve ", nColour, " from colour map ", cmap, " from color brewer, repeat ", nRep, " times\n"))
      
      ## Create colour vector:
      output <- rep(brewer.pal(nColour, cmap), nRep)
      output <- output[1:nLevel]
    } else {
      cat("Cannot find plotData, return colour red\n")
      output <- "red"
    }
  }
  
  cat(paste0("Retrieve colour scheme given variable name ", input, ": ", paste0(output, collapse = ", "), "\n"))
  return(output)
}

# ============================================================================ #
#### Add cue repetition: ####

add_cueRep <- function(data){
  #' Add new variables representing 
  #' (a) trial number per block;
  #' (b) cue repetition within block;
  #' (c) response on last trial with same cue.
  #' (d) outcome on last trial with same cue.
  #' @param data      data frame with trial-level data.
  #' @return data     same data frame with extra variables added.

  cat(paste0("Add trial number within block\n"))
  cat(paste0("Add cue repetition\n"))
  cat(paste0("Add response from last trial with same cue\n"))
  cat(paste0("Add outcome from last trial with same cue\n"))
  
  ## General settings:
  splitVar <- "block_n"
  stimVar <- "cue_block_n"
  respVar <- "response_n"
  outVar <- "outcome_n"
  
  ## Initialize variables:
  data$trialnr_block_n <- NA
  data$cueRep_n <- NA
  data$response_last_n <- NA
  data$outcome_last_n <- NA
  
  # --------------------------------------------------------------------- #
  ## Initialize counts:
  
  trialCount <- 1
  nCue <- 4
  cueCount <- rep(0, nCue)
  lastResp <- rep(NA, nCue)
  lastOut <- rep(NA, nCue)
  
  # --------------------------------------------------------------------- #
  ## First row:
  
  data$trialnr_block_n[1] <- trialCount
  
  thisCue <- data[1, stimVar] # identify cue
  cueCount[thisCue] <- cueCount[thisCue] + 1 # update
  data$cueRep_n[1] <- cueCount[thisCue] # save
  
  data$response_last_n[1] <- lastResp[thisCue] # save
  lastResp[thisCue] <- data[1, respVar] # update
  
  data$outcomelast_n[1] <- lastOut[thisCue] # save
  lastOut[thisCue] <- data[1, outVar] # update
  
  # --------------------------------------------------------------------- #
  ## All other rows:
  for (iRow in 2:nrow(data)){
    
    ## If different block: reset
    if(data[iRow, splitVar] != data[iRow - 1, splitVar]){
      trialCount <- 0
      cueCount <- rep(0, nCue)
      lastResp <- rep(NA, nCue)
    }
    
    ## Increment trial:
    trialCount <- trialCount + 1 # increment
    
    ## Increment cue repetition:
    thisCue <- data[iRow, stimVar] # identify cue
    cueCount[thisCue] <- cueCount[thisCue] + 1 # update
    
    ## Save:
    data$trialnr_block_n[iRow] <- trialCount # save
    data$cueRep_n[iRow] <- cueCount[thisCue] # save
    data$response_last_n[iRow] <- lastResp[thisCue] # save
    data$outcome_last_n[iRow] <- lastOut[thisCue] # save
    
    ## Update response and outcome:
    lastResp[thisCue] <- data[iRow, respVar] # update
    lastOut[thisCue] <- data[iRow, outVar] # update
    
  }  # end iRow 
  
  return(data)
}

# data[1:20, c("subject_n", "trialnr_n", "cue_n", "cueRep_n")]
# data[data$subject_n == 1 & data$cue_n == 2, c("subject_n", "trialnr_n", "cue_n", "cueRep_n", "response_n", "response_last_n")]
# data[data$subject_n == 1 & data$cue_n == 1, c("subject_n", "trialnr_n", "cue_n", "cueRep_n", "outcome_n", "outcome_last_n")]

# ============================================================================ #
#### Add lag 1 version of selected variable: ####

add_lag <- function(data, inputVar, nLag = 1){
  #' Overwrite reward lag 1 variable based on reward variable, delete after 
  #' choices and start of new round.
  #' @param data      data frame with trial-level data.
  #' @param inputVar  scalar string, name of variable to create lag for.
  #' @param nLag      scalar integer, lag to use (default: 1).
  #' @return data     same data frame with reward lag 1 variable corrected.
  
  # inputVar <- "outcome_n"; nLag = 1
  
  ## Variable to set new start of block:
  trialVar <- "trialnr_block_n"
  
  ## Create new variable name:
  outputVar <- inputVar
  outputVar <- gsub("_n", "", outputVar) # delete any trailing _n
  outputVar <- gsub("_f", "", outputVar) # delete any trailing _f
  outputVar <- gsub(" ", "", outputVar) # delete any spaces
  outputVar <- paste0(outputVar, "_lag", nLag)
  if (grepl("_n", inputVar, fixed = TRUE)){outputVar <- paste0(outputVar, "_n")}
  if (grepl("_f", inputVar, fixed = TRUE)){outputVar <- paste0(outputVar, "_f")}
  cat(paste0("Create new variable ", outputVar, "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Copy over from reward variable:
  
  data[, outputVar] <- NA # create empty variable
  data[(1 + nLag):nrow(data), outputVar] <- data[1:(nrow(data) - nLag), inputVar] # copy over with lag
  # data[1:20, c(inputVar, outputVar)]
  
  # -------------------------------------------------------------------------- #
  ### Delete for first nLag trials of new block:
  
  data[data[, trialVar] <= nLag, outputVar] <- NA # delete after leave choices
  
  # -------------------------------------------------------------------------- #
  ### Convert back to factor:
  
  if (grepl("_f", inputVar, fixed = TRUE)){data[, outputVar] <- factor(data[, outputVar])}
  
  # -------------------------------------------------------------------------- #
  ### Final NA assessment:
  
  cat(paste0("New variable ", outputVar, " has ", 100*round(mean(is.na(data[, outputVar])), 2), "% NAs\n"))
  
  return(data)
  
}

# ============================================================================ #
#### Recode character variables to proper factors ####

recode_char2fac <- function(data){
  #' Turn string variables into proper factors again
  #' @param data    data frame with trial-level data
  #' @return data   same data frame with all character variables turned into proper factors
  
  cat("Recode any character variable to factor\n")
  
  for (iCol in 1:ncol(data)){
    if (is.character(data[, iCol])){ # if character
      cat(paste0("Recode variable ", names(data)[iCol], "\n"))
      data[, iCol] <- factor(data[, iCol]) # recode to proper factor
    }
  }
  
  cat("Finished :-)\n")
  return(data)
}

# ============================================================================ #
#### Recode all numerical variables to factors ####

recode_num2fac <- function(data,variables=NULL){
  #' Recode selected numerical variable to factor
  #' @param data    data frame with trial-level data
  #' @variables     vector of strings, numerical variables to turn into factors (if not provided: all numerical variables in data frame)
  #' @return data   same data frame with all numerical variables also available as factor ("_f")
  
  cat("Recode selected numerical variables to factors\n")
  
  ### Determine if input variables given, otherwise take all numerical variables: 
  if(is.null(variables)){
    cat("No variables given, take all numerical variables\n")
    variables <- c()
    for (iCol in 1:ncol(data)){
      if (is.numeric(data[, iCol])){ # if numeric
        variables <- c(variables, names(data)[iCol]) # add variable name
      }
    }
  }
  
  ### Number variables:
  nVar <- length(variables)
  
  ### Loop through variables:  
  for (iVar in 1:nVar){ # loop through variables
    varName <- variables[iVar]
    cat(paste0("Recode variable ",varName, "\n"))
    newVarName <- paste0(varName, "_f") # new variable name
    data[, newVarName] <- factor(data[,varName]) # recode to proper factor
  }
  cat("Finished :-)\n")
  return(data)
}

# ============================================================================ #
#### Recode formula in Wilkinson notation to handle that can be used for saving models: ####

formula2handle <- function(formula){
  #' Create name based on formula with underscores instead of spaces without random-effects part to be used in plot names.
  #' @param formula   string, formula in Wilkinson notation.
  #' @return handle   string, converted to be without spaces but with underscores, with random-effects part removed.
  require(stringr)
  
  cat(paste0("Input: ", formula, "\n"))
  # https://stackoverflow.com/questions/38291794/extract-string-before
  # https://statisticsglobe.com/extract-substring-before-or-after-pattern-in-r
  
  handle <- formula # copy over
  
  # ------------------------------------------------------------------------ #
  ## For survival models: delete Surv() surrounding:
  
  if (grepl("Surv\\(", handle)){
    
    handle <- sub(".*Surv\\(", "", handle) # delete everything up until Surv(
    handle <- sub("\\)", "", handle) # delete first )
  }
  
  # ------------------------------------------------------------------------ #
  ## Extract until random effects parentheses:
  
  handle  <- sub("\\(.*", "", handle) # delete everything after (
  
  # ------------------------------------------------------------------------ #
  ## Delete only very last (!) plus before random-effects part:
  # https://stackoverflow.com/questions/44687333/remove-characters-after-the-last-occurrence-of-a-specific-character
  # handle <- sub("+[^_]+$", "", handle)
  
  if(grepl( "+", str_sub(handle, -2, -1), fixed = T)){
    handle <- substring(handle, 1, nchar(handle) - 3)
    
  }
  
  ## Replace every * by x:
  handle <- gsub("\\*", "x", handle) # extract until last plus
  
  ## Substitute every space with underscore:
  handle <- gsub(" ", "_", handle)
  
  cat(paste0("Output: ", handle, "\n"))
  return(handle)
  
}

# ============================================================================ #
#### Retrieve or fit & save model based on formula: ####

fit_lmem <- function(formula, useLRT = FALSE){
  #' Retrieve previously fitted model or fit model anew and save it.
  #' @param formula   string, formula in Wilkinson notation.
  #' @param useLRT    Boolean, compute p-values with LRTs using afex (TRUE) or not (FALSE; default).
  #' @return mod    fitted model object.

  require(lme4)
  require(afex)
  require(car)
  
  # ------------------------------------------------------------------------ #
  if (!exists("modData")){"modData does not exist"}

  # ------------------------------------------------------------------------ #
  ## Determine type of fitting:
  
  if (useLRT){
    fitType <- "LRT"
  } else {
    fitType <- "lme4"
  }
  
  # ------------------------------------------------------------------------ #
  ## Determine model family:
  
  DV <- sub("\\~.*", "", formula)
  if (grepl( "response", DV) | grepl( "ACC", DV) | grepl( "repeat", DV)){
    modFamily <- "binom"
  } else {
    modFamily <- "lin"
  } 
  
  ## Print specifics to console:
  cat(paste0("Fit model ", formula, " of family ", modFamily, " using ", fitType, "\n"))
  
  # ------------------------------------------------------------------------ #
  ## Determine model name:
  
  modName <- paste0(fitType, "_", modFamily, "_", formula2handle(formula), ".rds")
  
  # ------------------------------------------------------------------------ #
  ## Check if already exists; if yes, retrieve; if not, fit anew:
  
  if (file.exists(paste0(dirs$modelDir, modName))){
    
    cat(paste0(">>> Found ", modName, ", load \n"))
    mod <- readRDS(paste0(dirs$modelDir, modName))
    if (useLRT){
      print(anova(mod))
    } else {
      print(summary(mod))
      print(Anova(mod, type = "3")) # 1 p-value per factor
    }
    
  } else {
    
    # ---------------------------------------------------------------------- #
    ### Fit model:
    
    ## Start time:
    start.time <- Sys.time();
    cat(paste0(">>> Fit model with formula ", formula, "\n"))
    
    if (modFamily  == "binom"){ # if logistic regression
      if (useLRT){ # if LRT
        mod <- mixed(formula = formula, data = modData, method = "LRT", type = "III", family = binomial(), # all_fit = T,
                         control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
      } else { # if lme3
        mod <- glmer(formula = formula, data = modData, family = binomial(),
                     control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
      }
    } else { # if linear regression
      if (useLRT){ # if LRT
        mod <- mixed(formula = formula, data = modData, method = "LRT", type = "III", # all_fit = T,
                     control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
      } else { # if lme4
        mod <- lmer(formula = formula, data = modData,
                  control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
      }      
    }
    
    ## Stop time:
    end.time <- Sys.time(); beep()
    dif <- difftime(end.time,start.time); dif

    # ------------------------------------------------------------------------ #
    ### Print output:
    
    if (useLRT){
      print(anova(mod))
    } else {
      print(summary(mod))
      print(Anova(mod, type = "3")) # 1 p-value per factor
    }

    # ------------------------------------------------------------------------ #
    ### Save model:
    
    cat(paste0(">>> Save ", modName, "\n"))
    saveRDS(mod, paste0(dirs$modelDir, modName))
    cat(">>> Saved :-)\n")
    
  }
  
  return(mod)
  
}

# ============================================================================ #
#### Process variables (demean/ standardize/ conditional probabilities) for plotting in wide format: ####

prepare_data_plot_wide <- function(data, variables, centerSub = F, scaleSub = F, condProb = F){
  #' Extract variables, reshape, prepare such that they can easily be fed into raincloud plot
  #' @param data data frame, aggregated per subject, with variables \code{variables}
  #' @param variables vector of variable names which to plot
  #' @param centerSub whether to subtract overall mean per subject
  #' @param centerSub whether to z-standaridize data row per subject
  #' @param condProb whether to divide adjacent rows (1&2, 3&4, ...) by their sum to create conditional probabilities
  #' @return preprocessed data frame containing $x, $y, $xj, $subject, $subject_f
  
  # Select variables:
  subVar <- "subject_n"
  plotdata_wide <- data[ ,c(subVar, variables)]
  nVar <- length(variables)
  nRow <- nrow(plotdata_wide)
  nCol <- ncol(plotdata_wide)
  
  cat("Mean per variable:\n")
  print(colMeans(plotdata_wide[, 2:nCol], na.rm = T))
  
  if (centerSub == T){
    cat("Demean per subject\n")
    for (iRow in 1:nRow){ # iRow <- 1
      plotdata_wide[iRow, 2:nCol] <- plotdata_wide[iRow, 2:nCol] - mean(as.numeric(plotdata_wide[iRow, 2:nCol]), na.rm = T)
    } 
    # rowMeans(plotdata_wide[, 2:nCol]) # should all be zero
  }
  
  if (scaleSub == T){
    cat("Z-standardize per subject\n")
    for (iRow in 1:nRow){ # iRow <- 1
      plotdata_wide[iRow, 2:nCol] <- (plotdata_wide[iRow, 2:nCol] - mean(as.numeric(plotdata_wide[iRow, 2:nCol]), na.rm = T)) / std(as.numeric(plotdata_wide[iRow, 2:nCol]), na.rm = T)
    } 
  }
  
  # Conditional probabilities per condition:
  if (condProb == T){
    cat("Compute conditional probabilities for pairs of adjacent input variables\n")
    
    nPair <- nVar/2 # number adjacent pairs:
    for (iPair in 1:nPair){
      countVec <- plotdata_wide[, 1+2*iPair-1] + plotdata_wide[, 1+2*iPair] # add adjacent rows together
      plotdata_wide[, 1+2*iPair-1] <- plotdata_wide[, 1+2*iPair-1] / countVec # divide by total count
      plotdata_wide[, 1+2*iPair] <- plotdata_wide[, 1+2*iPair] / countVec # divide by total count
    }
    
    cat("New mean per variable:\n")
    print(colMeans(plotdata_wide[, 2:nCol], na.rm = T))
  }
  
  # Reshape to long format:
  plotdata_long <- reshape(plotdata_wide, idvar = subVar, varying = variables, 
                           timevar = "condition", v.names = "fixation", direction = "long")
  
  d <- plotdata_long
  
  # Make subject factor (e.g. when split per subject in plotting):
  d$subject_f <- factor(d$subject)
  
  # Assign to simple labels x and y:
  d$y <- d$fixation
  d$x <- d$condition
  
  # Add jitter to x (position):
  set.seed(321)
  d$j <- jitter(rep(0, nrow(d)), amount=.09)
  d$xj <- d$x + d$j
  
  # Return:
  return(d)
}

# ==================================================================================================== #
#### Process variables (demean/ standardize/ conditional probabilities) for plotting in long format ####

prepare_data_plot_long <- function(data, x, y, subVar = "subject_n", jitterNum = 0.09){
  require(plyr)
  
  if (!(x %in% names(data))){stop("Variable x not contained in data set")}
  if (!(y %in% names(data))){stop("Variable y not contained in data set")}
  
  data$x <- data[,x]
  data$y <- data[,y]
  data$subject <- data[,subVar]
  
  # data <- data[!(is.na(data$x)),] # delete data where IV is NA
  
  d <- ddply(data, .(x, subject), function(iData){
    
    y <- mean(iData$y, na.rm = T)
    
    return(data.frame(y))
    dev.off()})
  
  # Format x:
  d$condition <- d$x # save condition labels 
  d$x <- as.numeric(d$x) # turn numeric variable
  d$subject_f <- factor(d$subject)   
  
  # Add jitter to x (position):
  set.seed(321)
  if (jitterNum==0){ # if no jitter
    d$j <- 0
  } else { # if jitter
    d$j <- jitter(rep(0, nrow(d)), amount = jitterNum) # default: 0.9
  }
  d$xj <- d$x + d$j
  
  cat(paste0("Min = ",round(min(d$y, na.rm = T), 3), "; max = ",round(max(d$y, na.rm = T), 3), "\n"))
  
  return(d)
}

# ===============================================================================================#
#### Determine y-axis limits ####

determine_ylim_data_y <- function(data){
  #' Determine optimal y-axis limits based on some input heuristics
  #' @param data data frame, aggregated per subject, long format, with variable \code{y}
  #' @return yLim vector with to elements: minimal and maximal y-axis limit
  
  require(plyr)
  
  # Determine minimum and maximum:
  yMin <- min(data$y, na.rm = T)
  yMax <- max(data$y, na.rm = T)
  nRound <- 3
  cat(paste0("Detected yMin = ", round(yMin, nRound), ", yMax = ", round(yMax, nRound), "\n"))
  
  ## If binary input data:
  if(all(data$y[!is.nan(data$y)] %in% c(0, 1))){ 
    cat("Binary data, set y-axis limits to [0, 1]\n")
    yLim <- c(0, 1)
    
    ## if likely probability
  } else if (yMin >= 0 & yMax <= 1){ 
    cat("Set automatic yLim to [0, 1]\n")
    yLim <- c(0, 1)
    
    ## If positive number, but not huge:
  } else if(yMin >= 0 & yMax <= 10){ # if rather small positive number: plot from 0 onwards to 110% of yMax
    cat("Set automatic yLim to [0, round(yMax*1.1, 1)]\n")
    yLim <- c(0, round(yMax*1.1, 1))
    
    ## If positive number and huge:
  } else if ((yMin > 0) & (yMax > 100)) { # if very big number: plot from 0 onwards to next hundred
    cat("Set automatic yLim to [0, hundres]\n")
    yMin <- 0 # from zero onwards
    yMax <- round_any(yMax, 100, f = ceiling) # round up to next hundred
    yLim <- c(yMin, yMax)
    
    ## Else (if yMin negative): enforce symmetric yLim using the bigger of yMin and yMax
  } else { # take the numerically bigger one, symmetric
    cat("Set automatic yLim to be symmetric, use bigger of yMinAbs and yMaxAbs\n")
    yMaxAbs <- ceiling(c(abs(yMin), yMax)) # round to next integer
    yMaxAbs <- yMaxAbs[1] # only first entry
    yLim <- c(-1*yMaxAbs, 1*yMaxAbs)
  }
  
  ## Check if only 2 output elements:
  if (length(yLim) < 2){stop("yLim output has less than 2 elements, please check input\n")}
  if (length(yLim) > 2){stop("yLim output has more than 2 elements, please check input\n")}
  
  return(yLim)
} 

# ============================================================================ #
#### Plot single bar plot with individual data points: ####

custom_singlebar <- function(data, selVar, yLim = c(0, 1), isViolin = FALSE, hLine = NULL,
                             selCol = "grey80", xLab = "x", yLab = "y", main = NULL){
  #' Plot single bar or half-density (violin) with single points per data point.
  #' @param data data frame, with variable \code{selVar} to plot.
  #' @param selVar string, name of variable to plot.
  #' @param isViolin Boolean, plot half-density (violin) instead of bar (default: FALSE)
  #' @param selCol scalar string, color to use for violin/ bar (default: grey80).
  #' @param xLab string, label for x-axis (default: "x")
  #' @param yLab string, label for y-axis (default: "y")
  #' @param main string, title of plot (default: "NULL")
  #' @return prints and returns ggplot object.
  
  ## Load packages:
  require(ggplot2)
  require(gghalves)
  require(ggbeeswarm)
  
  ## Aggregate again with Rmisc:
  # library(Rmisc)
  # summary_d <- summarySEwithin(d,measurevar = "ACC", idvar = "subject_n", na.rm = T)
  
  if(!(selVar %in% names(data))){stop(paste0("Variable ", selVar, " not found in data"))}
  
  # -------------------------------------------------------------------------- #
  ### Fixed settings:
  
  FTS <- retrieve_plot_defaults("FTS") # 20 # 30
  colAlpha <- .70
  LWD <- retrieve_plot_defaults("LWD")
  
  # -------------------------------------------------------------------------- #
  ### Prepare data set:
  
  ## Add jittered x-axis position:
  nSub <- nrow(data)
  data$x <- rep(1, nSub)
  data$j <- jitter(rep(0, nrow(data)), amount = .05)
  data$xj <- data$x + data$j
  
  ## Repeat selected variable:
  data$y <- data[, selVar]
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  p <- ggplot(data = data, aes(x = x, y = y)) # initialize
  
  # -------------------------------------------------------------------------- #
  ## Bar or violin:
  if (isViolin){ # if violin
    
    p <- p + geom_half_violin(color = "black", fill = selCol, alpha = colAlpha, trim = FALSE)
    
  } else { # if bar
    
    p <- p + stat_summary(fun = mean, geom = "bar", fill = selCol, alpha = colAlpha,
                          color = 'black', width = 0.3, lwd = LWD)
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Confidence intervals: 
  
  # p <- p + stat_summary(fun.data = mean_cl_normal, geom =
  #                         "errorbar", width = 0.05, lwd = LWD)
  
  # -------------------------------------------------------------------------- #
  ## Standard error: 
  
  ## Compute data summaries:
  summary_d <- data.frame(x = 1)
  summary_d$y <- mean(data$y, na.rm = T)
  summary_d$sd <- sd(data$y, na.rm = T)
  summary_d$se <- se(data$y, na.rm = T)
  SEweight <- 1
  
  ## Add error bar manually:
  p <- p + geom_errorbar(data = summary_d, aes(y = y, ymin = y - se * SEweight, ymax = y + se * SEweight), # x = x, 
                         width = 0.05, lwd = LWD) # , color = "black", alpha = 1)
  
  # -------------------------------------------------------------------------- #
  ## Individual data points:
  # p <- p + geom_beeswarm(color = "black", fill = color, alpha = colAlpha)
  p <- p + geom_point(data = data, aes(x = xj), color = "black", fill = "grey40",
                      shape = 21, size = 4,
                      alpha = colAlpha)
  
  # -------------------------------------------------------------------------- #
  ## Other settings:
  
  if (!(is.null(hLine))){ # if conditional probabilities:
    p <- p + geom_hline(yintercept = hLine, linetype = 2, color = "black")
  }
  # if (mean(yLim) == 0.5){ # if conditional probabilities:
  #   p <- p + geom_hline(yintercept=0.5, linetype=2, color = "black")
  # }
  
  ## X-axis:
  p <- p + scale_x_continuous(limits = c(0.5, 1.5), breaks = c(0, 1, 2), labels = c("", "", ""))
  
  # Y-axis:
  if (yLim[1] >= 0 & yLim[-1] <= 1){
    p <- p + scale_y_continuous(limits = yLim, breaks = seq(0, 1, 0.25)) # steps of 0.25
  } else if (yLim[1] >= 0 & yLim[-1] <= 5){
    p <- p + scale_y_continuous(limits = yLim, breaks = seq(0, ceiling(yLim[-1]), 0.50)) # steps of 0.50
  } else if (yLim[1] >= 0 & yLim[-1] <= 10){
    p <- p + scale_y_continuous(limits = yLim, breaks = seq(0, ceiling(yLim[-1]), 1)) # steps of 0.50
  } else {
    # p <- p + coord_cartesian(ylim = yLim)
    p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1], yLim[-1], 5)) # steps of 1
  }
  
  # Axis labels:
  p <- p + xlab(xLab) + ylab(yLab)
  
  ## Add title:
  if(!is.null(main)){p <- p + ggtitle(main)}    
  
  ## Add theme:
  p <- p + theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5), # center title 
                 legend.text = element_text(size = FTS))
  
  # Print plot in the end:
  print(p)
  return(p)
  
}

# =============================================================================================== #
#### REGRESSION LINES 1 IV: Plot regression line per group and per subject based on model output: #####

custom_regressionline1 <- function(mod, selEff, xLim = NULL, yLim = NULL, useEffect = TRUE, xVec = NULL,
                                   selCol = "red", margin = NULL, xLab = "x", yLab = "y", main = NULL, FTS = NULL){
  #' Plot group-level regression line and subject-level regression lines based on 1 continuous predictor.
  #' @param mod model fitted with lme4
  #' @param selEff string, name of predictor to plot
  #' @param xLim vector of two numbers for y-axis limits
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data)
  #' @param subVar string, name of grouping variable (default: subject)
  #' @param useEffect boolean, extract upper and lower bound of confidence interval from effects(mod) (TRUE) or compute yourself based on vcov(mod) (FALSE)
  #' @param xVec vector of two numbers for x-axis ticks, to be used only if useEffect = FALSE.
  #' @param selCol strings (HEX colors), colors for line and error shade (default: red)
  #' @param margin vector of 4 numbers, margin of plot (default: NULL)
  #' @param xLab string, label for x-axis (default: "x")
  #' @param yLab string, label for y-axis (default: "y")
  #' @param main string, title of plot (default: "Plot")
  #' @param FTS integer, font size for axes ticks and labels and title.
  #' @return makes regression line plot
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.Rmd
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.pdf
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  require(ggplot2)
  require(lme4)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## General settings:
  
  colAlpha <- .95 # transparency
  LWD <- retrieve_plot_defaults("LWD") # line width
  
  if (is.null(FTS)){
    ## Font sizes for ordinary viewing: 15
    # FTS <- 15
    ## Font sizes for saving: 30
    FTS <- retrieve_plot_defaults("FTS")
    cat(paste0("No font size provided, use font size ", FTS), "\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  
  ## Extract group-level and subject-level data sets:
  groupCoefs <- fixef(mod) # fixed effects
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects (only random effects for first structure)
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  if (sum(grepl(selEff, effNames, fixed = T)) == 0){stop("selEff not a predictor in mod")}
  selEff1 <- effNames[grep(selEff, effNames, fixed = T)]
  selEff1 <- selEff1[1] # first effect
  
  ## Check presence of interactions:
  if(sum(grepl(":", effNames, fixed = T)) > 0){warning("Model contains interactions; predictions might be off, especially if size of interaction effect is substantial")}
  
  ## Locate effect, count subjects:
  iCol <- which(names(subCoefs) == selEff1) # localize where in subCoefs effect of interest is
  allCols <- 1:length(subCoefs)
  otherCols <- allCols[allCols != iCol]
  nSub <- nrow(subCoefs)
  
  # --------------------------------------------------------------------- #
  ## Generate x-axis values for which to generate plots:
  
  tmp <- effect(selEff, mod) # retrieve objects from effects
  
  if (useEffect){ # if x samples automatically generated by effects() to be used
    
    xVecEval <- as.numeric(tmp$model.matrix[, iCol])
    
  } else { # if use by hand
    
    if (is.null(xVec)){ # if no xVec provided
      
      cat("No x-axis samples provided, extract from effects(mod)\n")
      xVecEval <- as.numeric(tmp$model.matrix[, iCol])
      
    } else {
      
      xVecEval <- seq(xVec[1], xVec[2], (xVec[2] - xVec[1]) / 100) # 100 samples between min and max
      
    }
    
  }
  
  xLen <- length(xVecEval) # number of x-axis samples
  
  if (is.null(xLim)){ # if no x-axis limits provided
    xLim <- c(min(xVecEval), max(xVecEval))
  }
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot; also assign xVecEval to y
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines + ribbons:
  cat("Draw random-effects lines\n")
  
  inputData <- mod@frame # extract input data
  subIdx <- inputData[, length(inputData)] # extract subject indices to compute marginal means
  subVec <- unique(subIdx)
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    subID <- subVec[iSub]
    margMean <- colMeans(model.matrix(mod)[subIdx == subID, ]) # mean of each regressor in design matrix
    iIntercept <- as.numeric(as.numeric(subCoefs[iSub, otherCols]) %*% as.numeric(margMean[otherCols])) # intercept is estimate of y at mean of all other regressors (marginal effect)
    # iIntercept <- subCoefs[iSub, 1] # extract intercept --> works only if other variables centered or very small
    iSlope <- subCoefs[iSub, iCol] # extract slope
    yVecEval <- iIntercept + xVecEval * iSlope
    if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # bring to response scale
    
    ## Create single data frame (2 points should be enough to draw line, i.e. xmin and xmax):
    d <- data.frame(x = xVecEval, y = yVecEval)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + geom_path(data = d, aes(x = x, y = y), color = 'grey40', # color = 'grey70'
                       alpha = 0.35, linewidth = 1.0)
    # Used till 2021-11-21: grey50, size = 1
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level line:
  cat("Draw fixed-effects line\n")
  
  ## Compute group-level predicted variables:
  
  # iIntercept <- as.numeric(groupCoefs[1]) # extract intercept
  margMean <- colMeans(model.matrix(mod)) # mean of each regressor in design matrix
  iIntercept <- as.numeric(groupCoefs[otherCols] %*% margMean[otherCols]) # intercept is estimate of y at mean of all other regressors (marginal effect)
  # iIntercept <- groupCoefs[1] # extract intercept --> works only if other variables centered or very small
  iSlope <- as.numeric(groupCoefs[iCol]) # extract slope
  yVecEval <- iIntercept + xVecEval * iSlope # recompute before transform
  if (isGLMM(mod)){yVecEval <- mod@resp$family$linkinv(yVecEval)} # bring to response scale
  
  ## Based on effects output:
  # yVecEval2 <- as.numeric(tmp$fit) # y-axis coordinates (untransformed)
  # if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # bring to response scale
  
  
  ## Create single data frame (2 points should be enough to draw line, i.e. xmin and xmax):
  d <- data.frame(x = xVecEval, y = yVecEval)
  
  ## Thick line connecting means (plot line first and points on top):
  p <- p + geom_path(data = d, aes(x = x, y = y), color = selCol, linewidth = LWD)
  
  # ----------------------------------------------------- #
  ## Error shades:
  if (useEffect){
    
    # ------------------------------- #
    ## Option A: Extract from effect() object:
    
    tmp <- effect(selEff, mod)
    
    # Mean/ line itself:
    yVecEval <- as.numeric(tmp$fit) # y-axis coordinates (untransformed)
    if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # transform to response scale
    
    # Lower and upper limit of CI interval:
    ymin <- t(tmp$lower) # lower bound of CI interval
    ymax <- t(tmp$upper) # upper bound of CI interval
    if (isGLMM(mod)){ymin <-mod@resp$family$linkinv(ymin)} # transform to response scale
    if (isGLMM(mod)){ymax <-mod@resp$family$linkinv(ymax)} # transform to response scale
    
  } else {
    
    # ------------------------------- #
    ## Option B: Compute SE yourself:
    # https://github.com/cran/effects/blob/master/R/Effect.R
    
    mmat <- matrix(0, xLen, ncol(subCoefs)) # initialize design matrix: by default all regressors 0
    mmat[, 1] <- rep(1, xLen) # add intercept
    mmat[, iCol] <- xVecEval # add regressor of interest
    V <- vcov(mod, complete=FALSE) # covariance matrix from model
    vcov <- mmat %*% V %*% t(mmat) # multiply with design matrix
    var <- diag(vcov) # variance
    se <- sqrt(var) # se # see tmp$se
    conf <- 1.96
    
    # iIntercept and iSlope computed above:
    yVecEval <- iIntercept + xVecEval * iSlope # recompute before transform
    ymin <- yVecEval - se*conf # lower bound of CI interval
    ymax <- yVecEval + se*conf # upper bound of CI interval
    if (isGLMM(mod)){ymin <-mod@resp$family$linkinv(ymin)} # bring to response scale
    if (isGLMM(mod)){ymax <-mod@resp$family$linkinv(ymax)} # bring to response scale
    # ymin;ymax
  }
  
  # ------------------------------- #
  ## Plot error bars/ shades:
  
  d <- data.frame(x = xVecEval, y = yVecEval, ymin = ymin, ymax = ymax)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol, alpha = 0.20, linetype = 0)
  
  # --------------------------------------------------------------------- #
  cat("Adjust axes, labels\n")
  
  # Y-axis:
  p <- p + coord_cartesian(xlim = xLim) 
  if (!is.null(yLim)){p <- p + coord_cartesian(ylim = yLim)} 
  
  # X-axis:
  xTickVec <- round(seq(xLim[1], xLim[2],(xLim[2] - xLim[1]) / 2), 2)
  p <- p + scale_x_continuous(breaks=xTickVec, labels=xTickVec)
  
  # Labels:
  p <- p + xlab(xLab) + ylab(yLab)
  
  # Other settings:
  if (!is.null(main)){
    cat("Add title\n")
    p <- p + ggtitle(main)
  }
  
  p <- p + theme_classic()
  if (!is.null(margin)){
    cat("Adjust margin\n")
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5), 
                 legend.text = element_text(size = FTS))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### REGRESSION LINES 2 IV: Plot regression line per condition per group and per subject based on model output: #####

custom_regressionline2 <- function(mod, xVar, zVar, 
                                   xLim = NULL, yLim = NULL, xVec = NULL,
                                   selCol = NULL, margin = NULL, 
                                   xLab = NULL, yLab = NULL, main = NULL, FTS = NULL){
  #' Plot group-level regression line and subject-level regression lines based on 1 continuous predictor and 1 binary predictor.
  #' @param mod model fitted with lme4.
  #' @param xVar string, name of continuous predictor to plot on x-axis.
  #' @param zVar string, name of binary predictor to plot with different colors.
  #' @param subVar string, name of grouping variable (default: subject).
  #' @param xLim vector of two numbers for y-axis limits.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @param xVec vector of two numbers for x-axis ticks.
  #' @param selCol vector of strings (HEX colors), colors for line and error shade (default: retrieved from retrieve_colour()).
  #' @param margin vector of 4 numbers, margin of plot (default: NULL).
  #' @param xLab string, label for x-axis (default: retrieved from substitute_label()).
  #' @param yLab string, label for y-axis (default: retrieved from substitute_label()).
  #' @param main string, title of plot (optional).
  #' @param FTS integer, font size for axes ticks and labels and title.
  #' @return makes regression line plot.
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.Rmd
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.pdf
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  
  require(ggplot2)
  require(lme4)
  
  # --------------------------------------------------------------------- #
  ## General settings:
  
  alphaSub <- 0.10
  alphaShade <- 0.20
  sizeGroup <- 1.5
  sizeSub <- 1
  LWD <- retrieve_plot_defaults("LWD")
  
  if (is.null(FTS)){
    ## Font sizes for ordinary viewing: 15
    # FTS <- 15
    ## Font sizes for saving: 30
    FTS <- retrieve_plot_defaults("FTS")
    cat(paste0("No font size provided, use font size ", FTS), "\n")
  }
  
  if(is.null(selCol)){selCol <- retrieve_colour(zVar)}
  if(is.null(xLab)){xLab <- substitute_label(xVar)}
  if(is.null(yLab)){yLab <- substitute_label(all.vars(terms(mod))[1])}
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  
  ## Extract group-level and subject-level data sets:
  groupCoefs <- fixef(mod) # fixed effects
  
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects (only random effects for first structure)
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  
  cat("Extract predictor names from model\n")
  if (sum(grepl(xVar, effNames, fixed = T)) == 0){stop("xVar not a predictor in mod")}
  xVar1 <- effNames[grep(xVar, effNames, fixed = T)] # check version in effect names
  xVar1 <- xVar1[1] # only first match
  
  if (sum(grepl(zVar, effNames, fixed = T)) == 0){stop("zVar not a predictor in mod")}
  zVar1 <- effNames[grep(zVar, effNames, fixed = T)] # check version in effect names
  zVar1 <- zVar1[1] # only first match
  
  isRev = FALSE
  intVar <- paste0(xVar, ":", zVar) # combine as name should look like (without 1)
  intVar1 <- paste0(xVar1, ":", zVar1) # combine as name should look like (with 1)
  cat(paste0("Assume interaction is ", intVar, "\n"))
  intVar1 <- effNames[grep(intVar1, effNames, fixed = T)] # check version in effect names (overwrite)
  intVar1 <- intVar1[1] # select first one
  if (is.na(intVar1)){ # if doesn't exist: try reverse
    isRev = TRUE
    intVar <- paste0(zVar, ":", xVar) # combine as name should look like (without 1)
    intVar1 <- paste0(zVar1, ":", xVar1) # combine as name should look like (with 1)
    cat(paste0("Assume interaction is ", intVar, "\n"))
    intVar1 <- effNames[grep(intVar1, effNames, fixed = T)] # check version in effect names (overwrite)
    intVar1 <- intVar1[1] # select first one
    if (is.na(intVar1)){stop("Interaction not a predictor in mod")} 
  }
  
  ## Indices of effects:
  xCol <- which(names(subCoefs) == xVar) # localize where in subCoefs effect of interest is
  zCol <- which(names(subCoefs) == zVar1) # localize where in subCoefs effect of interest is
  intCol <- which(names(subCoefs) == intVar1) # localize where in subCoefs effect of interest is
  # cat(paste0("xCol = ", xCol, ", zCol = ", zCol, ", intCol = ", intCol, "\n"))
  
  ## Number of subjects:
  nSub <- nrow(subCoefs)
  
  ## X-axis for which to generate plots:
  tmp <- effect(intVar, mod) # retrieve objects from effects (without the 1)
  if (is.null(xVec)){
    xVecEval <- sort(unique(as.numeric(tmp$model.matrix[, xCol])))
  } else {
    xVecEval <- xVec # copy over input
  }
  
  xLen <- length(xVecEval) # number of x-axis samples
  
  if (is.null(xLim)){ # if no x-axis limits provided
    xLim <- c(xVecEval[1], xVecEval[xLen])
  }
  
  ## Z-axis with how levels of factor are coded:
  zVecEval <- unique(as.numeric(tmp$model.matrix[, zCol]))
  if(length(zVecEval) > 2){stop("zVec is variable with > 2 levels, must have only 2 levels")}
  
  # --------------------------------------------------------------------- #
  ## Initialize empty ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot; assign xVecEval also to y
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines:
  
  cat("Draw random-effects lines\n")
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    ## Intercepts:
    iInter1 <- subCoefs[iSub, 1] + zVecEval[1] * subCoefs[iSub, zCol] # extract intercept condition 1
    iInter2 <- subCoefs[iSub, 1] + zVecEval[2] * subCoefs[iSub, zCol] # extract intercept condition 2
    
    ## Slopes:
    iSlope1 <- subCoefs[iSub, xCol] + zVecEval[1] * subCoefs[iSub, intCol] # extract slope condition 1
    iSlope2 <- subCoefs[iSub, xCol] + zVecEval[2] * subCoefs[iSub, intCol] # extract slope condition 2
    
    # Simulated y-data per subject:
    yVecEval1 <- iInter1 + xVecEval * iSlope1
    yVecEval2 <- iInter2 + xVecEval * iSlope2
    
    if (isGLMM(mod)){yVecEval1 <- mod@resp$family$linkinv(yVecEval1)} # bring to response scale
    if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # bring to response scale
    
    ## Create single data points (2 should be enough, i.e. xmin and xmax)
    d <- data.frame(x = xVecEval, y1 = yVecEval1, y2 = yVecEval2)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + 
      geom_path(data = d, aes(x = x, y = y1), color = selCol[1],
                alpha = alphaSub, linewidth = sizeSub) + 
      geom_path(data = d, aes(x = x, y = y2), color = selCol[2],
                alpha = alphaSub, linewidth = sizeSub)
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level line:
  
  cat("Draw fixed-effects line\n")
  
  ## Extract from effect() object:
  tmp <- effect(intVar, mod)
  
  ## Extract group-level line + CIs from model: 
  
  ## Create indices:
  if (isRev == TRUE){ # if continuous regressor comes first
    idx1 <-seq(1, xLen*2, 2) # odd indices
    idx2 <-seq(2, xLen*2, 2) # even indices
  } else { # if continuous regressor comes second
    idx1 <- 1:xLen # first half of predicted values
    idx2 <- (xLen+1):(2*xLen) # second half of predicted values
  }
  
  ## Line itself (means):
  yVecEval1 <- as.numeric(tmp$fit)[idx1] # y-axis coordinates (untransformed)
  yVecEval2 <- as.numeric(tmp$fit)[idx2] # y-axis coordinates (untransformed)
  
  ## Transform:
  if (isGLMM(mod)){yVecEval1 <- mod@resp$family$linkinv(yVecEval1)} # transform to response scale
  if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # transform to response scale
  
  ## Lower and upper limit of CI interval:
  ymin1 <- t(tmp$lower)[idx1] # lower bound of CI interval condition 1
  ymin2 <- t(tmp$lower)[idx2] # lower bound of CI interval condition 2
  ymax1 <- t(tmp$upper)[idx1] # upper bound of CI interval condition 1
  ymax2 <- t(tmp$upper)[idx2] # upper bound of CI interval condition 2
  
  ## Transform:
  if (isGLMM(mod)){ymin1 <- mod@resp$family$linkinv(ymin1)} # transform to response scale
  if (isGLMM(mod)){ymin2 <- mod@resp$family$linkinv(ymin2)} # transform to response scale
  if (isGLMM(mod)){ymax1 <- mod@resp$family$linkinv(ymax1)} # transform to response scale
  if (isGLMM(mod)){ymax2 <- mod@resp$family$linkinv(ymax2)} # transform to response scale
  
  # --------------------------------------------------------------------- #
  ## Thick line connecting group-level means:
  
  d <- data.frame(x = xVecEval, y1 = yVecEval1, y2 = yVecEval2)
  p <- p + 
    geom_path(data = d, aes(x = x, y = y1), color = selCol[1], linewidth = sizeGroup) + 
    geom_path(data = d, aes(x = x, y = y2), color = selCol[2], linewidth = sizeGroup)
  
  # --------------------------------------------------------------------- #
  ## Error shades:
  
  ## Plot error bars/ shades:
  d <- data.frame(x = xVecEval, y = yVecEval1, ymin = ymin1, ymax = ymax1)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol[1], alpha = alphaShade, linetype = 0)
  d <- data.frame(x = xVecEval, y = yVecEval2, ymin = ymin2, ymax = ymax2)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol[2], alpha = alphaShade, linetype = 0)
  
  # --------------------------------------------------------------------- #
  ## Further plot settings:
  
  cat("Adjust axes, labels\n")
  
  ## Y-axis:
  p <- p + coord_cartesian(xlim = xLim, ylim = yLim) 
  
  ## X-axis:
  # xTickVec <- round(seq(xLim[1], xLim[2], (xLim[2] - xLim[1]) / 2), 2) # extremes, middle
  xTickVec <- round(xVecEval, 2) # extremes, middle
  p <- p + scale_x_continuous(breaks=xTickVec, labels=xTickVec)
  
  ## Labels:
  p <- p + xlab(xLab) + ylab(yLab)
  
  ## Title:
  if (!is.null(main)){
    cat("Add title\n")
    p <- p + ggtitle(main)
  }
  
  ## Theme:  
  p <- p + theme_classic()
  
  ## Margin:
  if (!is.null(margin)){
    cat("Adjust margin\n")
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5), 
                 legend.text = element_text(size = FTS))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### REGRESSION BARS 1 IV: Plot regression bars per group and per subject based on model output: #####

custom_regressionbar1 <- function(mod, selEff, selCol = "red", 
                                  xLab = "x", yLab = "y", main = NULL, xLabels = NULL,
                                  margin = NULL, FTS = NULL, yLim = NULL){
  #' Plot group-level regression bar and subject-level regression lines based on 1 binary predictor.
  #' @param mod model fitted with lme4.
  #' @param selEff string, name of predictor to plot.
  #' @param subVar string, name of grouping variable (default: subject).
  #' @param selCol strings (HEX colors), colors for line and error shade (default: "red" for all).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (default: NULL).
  #' @param xLabels vector of strings, x-axis ticks (optional; otherwise retrieved from model).
  #' @param margin vector of 4 numbers, margin of plot (default: NULL).
  #' @param FTS integer, font size for axes ticks and labels and title.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @return makes regression line plot
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  
  require(ggplot2)
  require(lme4)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # --------------------------------------------------------------------- #
  ## General settings:
  
  colAlpha <- .95
  LWD <- retrieve_plot_defaults("LWD")
  
  if (is.null(FTS)){
    ## Font sizes for ordinary viewing: 15
    # FTS <- 15
    ## Font sizes for saving: 30
    FTS <- retrieve_plot_defaults("FTS")
    cat(paste0("No font size provided, use font size ", FTS), "\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  groupCoefs <- fixef(mod) # fixed effects
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  if (sum(grep(selEff, effNames, fixed = T)) == 0){stop("selEff not a predictor in mod")}
  selEff1 <- effNames[grep(selEff, effNames, fixed = T)] # check version in effect names
  selEff1 <- selEff1[1] # only first match
  
  ## Check presence of interactions:
  if(sum(grepl(":", effNames, fixed = T)) > 0){
    useEffect = FALSE
    warning("Model contains interactions; set useEffect to FALSE, check if predictions make sense")
  }
  
  ## Locate effect:
  iCol <- which(names(subCoefs) == selEff1) # localize where in subCoefs effect of interest is
  
  ## Count subjects:
  nSub <- nrow(subCoefs)
  
  ## x-coordinates for which to plot:
  tmp <- effect(selEff, mod) # retrieve objects from effects
  xVecEval <- as.numeric(tmp$model.matrix[, iCol])
  xLen <- length(xVecEval)
  
  ## In case of sum-to-zero coding: Flip everything
  sum2zero <- F
  if(all(xVecEval == c(1, -1))){sum2zero <- T; message("Sum-to-zero coding detected, flip y values, colors, x-axis labels")}
  
  ## Complete color vector:
  selCol <- rep(selCol, length.out = xLen)
  if (sum2zero){selCol <- rev(selCol)} # flip of sum-to-zero coding
  
  ## X-axis tick labels:
  if (is.null(xLabels)){
    xLabels <- tmp$variables[[1]]$levels
  }
  if (sum2zero){xLabels <- rev(xLabels)} # flip of sum-to-zero coding
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines + ribbons:
  cat("Draw random-effects lines\n")
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    iInter <- subCoefs[iSub, 1] # extract intercept
    iSlope <- subCoefs[iSub, iCol] # extract slope
    yVecEval <- iInter + xVecEval * iSlope # swap both x-axis positions
    if (isGLMM(mod)){yVecEval <- mod@resp$family$linkinv(yVecEval)} # bring to response scale
    if (sum2zero){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
    
    ## Create single data frame (2 data points should be enough, i.e. xMin and xMax):
    d <- data.frame(x = xVecEval, y = yVecEval)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + geom_path(data = d, aes(x = x, y = y), color = 'grey40', # color = 'grey70'
                       alpha = 0.35, linewidth = 1)
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level mean and error-bar:
  cat("Draw fixed-effects line\n")
  
  # groupSE <- VarCorr(mod)
  # iCol <- which(colnames(groupCoefs)==selEff) # localize where in groupCoefs effect of interest is
  iInter <- as.numeric(groupCoefs[1]) # extract intercept
  iSlope <- as.numeric(groupCoefs[iCol]) # extract slope
  yVecEval <- iInter + xVecEval * iSlope
  if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # bring to response scale
  if (sum2zero){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
  
  ## Create single data frame (2 data points should be enough, i.e. xMin and xMax):
  d <- data.frame(x = xVecEval, y = yVecEval)
  
  # ------------------------------- #
  ## Thick line connecting means (plot line first and points on top):
  p <- p + geom_path(data = d, aes(x = x, y = y), color = "black", linewidth = LWD)
  
  # ------------------------------- #
  ## Point for mean:
  p <- p + geom_point(data = d, aes(x = x, y = y), # point
                      color = selCol, alpha = 1, size = 5) # size = 2
  
  # ------------------------------- #
  ## Error shades:
  
  ymin <- tmp$lower 
  ymax <- tmp$upper 
  if (isGLMM(mod)){ymin <- mod@resp$family$linkinv(ymin)} # bring to response scale
  if (isGLMM(mod)){ymax <- mod@resp$family$linkinv(ymax)} # bring to response scale
  if (all(xVecEval == c(1, -1))){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
  if (all(xVecEval == c(1, -1))){ymin <- rev(ymin)} # flip of sum-to-zero coding
  if (all(xVecEval == c(1, -1))){ymax <- rev(ymax)} # flip of sum-to-zero coding
  d <- data.frame(x = xVecEval, y = yVecEval, ymin = ymin, ymax = ymax)
  p <- p + geom_errorbar(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                         color = selCol, width = 0.15, size = 1.5, alpha = .6)
  
  # --------------------------------------------------------------------- #
  ## Further plot settings:
  
  cat("Adjust axes, labels\n")
  
  ## Y-axis:
  if (!is.null(yLim)){
    p <- p + coord_cartesian(ylim = yLim) 
    if (yLim[1] == 0 & yLim[2] == 1){
      cat("Y-axis from 0 to 1\n")
      p <- p + scale_y_continuous(limits = yLim, breaks = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
    }
  }
  
  ## X-axis:
  p <- p + coord_cartesian(xlim = c(min(xVecEval) - 0.5, max(xVecEval) + 0.5)) 
  p <- p + scale_x_continuous(breaks = xVecEval, labels = xLabels)
  
  ## Labels:
  p <- p +  xlab(xLab) + ylab(yLab)
  
  ## Title:
  if (!is.null(main)){
    p <- p + ggtitle(main) # title off for printing for poster
  }
  
  ## Theme:
  p <- p + theme_classic() # theme
  
  ## Margin:
  if (!is.null(margin)){
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5), # center title 
                 legend.text = element_text(size = FTS))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### Plot group-level and subject-level coefficients as dots in horizontal dot-plot: #####

custom_coefplot <- function(mod, plotSub = TRUE, plotText = FALSE, dropIntercept = FALSE, revOrder = FALSE,
                            xLab = "Regression weight", yLab = "Predictor", main = NULL,
                            selCol = "blue", yLabels = NULL, xLim = NULL, FTS = NULL){
  #' Plot group-level  and subject-level coefficients as dots in horizontal dot-plot.
  #' @param mod model fitted with lme4.
  #' @param plotSub Boolean, whether to plot per-subject effects (TRUE) or not (FALSE; default: TRUE).
  #' @param plotText Boolean, whether to print value of group-level coefficient next to dot (TRUE) or not (FALSE; default: FALSE).
  #' @param dropIntercept Boolean, do not plot intercept (TRUE; default: false).
  #' @param revOrder Boolean, revert order of predictors (first one on top, TRUE) or not (last one on top FALSE) (default: false).
  #' @param xLab string, label for x-axis (default: "Regression weight").
  #' @param yLab string, label for y-axis (default: "Predictor").
  #' @param main string, title of plot (default: NULL).
  #' @param selCol strings (HEX colors), colors for group-level dots and error lines (default: "blue" for all).
  #' @param yLabels vector of strings, y-axis ticks (default: terms extracted from mod).
  #' @param xLim vector of two numerics, x-axis limits (optional).
  #' @param FTS scalar integer, font size (optional; default: NULL).
  #' @return coefplot created with ggplot.
  
  # -------------------------------------------------------------------------- #
  ## Required packages:
  
  require(ggplot2)
  require(lme4)
  require(arm) # for se.fixef
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Fixed settings:
  
  SEweight <- 1.96 # width of whiskers (1.96 for two-sided 95%-CI)
  groupDotSize <- 5 # size of fixed-effects dot; used to be 5
  subDisplacement <- -0.20 # systematic vertical displacement of per-subject dots
  subJitter <- 0.05 # amount of jitter added to per-subject dots; used to be 0.07
  textOffset <- 0.35 # vertical upwards displacement of text; used to be 0.3
  lineWidth <- retrieve_plot_defaults("LWD") # linewidth of axes
  if(is.null(FTS)){
    FTS <- retrieve_plot_defaults("FTS") # font size for all text: 30 or 15 
  }
  colAlpha <- 0.6 # transparency of per-subject dots
  nRound <- 3 # how much to round plotted text.
  
  # -------------------------------------------------------------------------- #
  ## Extract group-level information from input:
  
  ## a) If mixed effects model: 
  modClass <- class(mod)
  if(modClass %in% c("glmerMod", "lmerMod", "lmerTest", "lmerModLmerTest")){ # from lme4
    
    # Extract fixed effect:
    meanVec <- as.numeric(fixef(mod))
    seVec <- se.fixef(mod)
    
    if (is.null(yLabels)){ # if not provided
      labelsVec <- colnames(mod@pp$X)
      labelsVec <- substitute_label(labelsVec) # translate labels
    } else { # if provided
      labelsVec <- yLabels
    }
    
    ## Concatenate group-level values to data frame:
    groupCoefs <- data.frame(labelsVec, meanVec, seVec)
    names(groupCoefs) <- c("label", "mean", "se")
    
    ## Subject-level coefficients:
    subCoefs <- coef(mod)[[1]]
    nSub <- nrow(subCoefs)
    
    ## b) If flat regression:
  } else if (is.list(mod) & "bMat" %in% names(mod)) {
    
    # Compute mean and SD of per-subject coefficients:
    meanVec <- colMeans(mod$bMat, na.rm = T)
    nSub <- nrow(mod$bMat)
    seVec <- as.numeric(sapply(data.frame(mod$bMat), sd, na.rm = T))/sqrt(nSub)
    
    if (is.null(yLabels)){ # if not provided
      labelsVec <- names(coef(mod$modList[[1]]))
      labelsVec <- substitute_label(labelsVec) # translate labels
    } else { # if provided
      labelsVec <- yLabels
    }
    
    ## Concatenate group-level values to data frame:
    groupCoefs <- data.frame(labelsVec, meanVec, seVec)
    names(groupCoefs) <- c("label", "mean", "se")
    
    ## Subject-level coefficients:
    subCoefs <- data.frame(mod$bMat)
    names(subCoefs) <- names(coef(mod$modList[[1]])) # copy over regressor names
    nSub <- nrow(subCoefs)
    
  } else {
    stop("Unknown input")
  }
  
  ## Text labels with significance stars based on z-values:
  # (1 - pnorm(1.64))*2
  # (1 - pt(1.96, df = 1000))*2
  # nRound <- 3
  ## Compute absolute z-value for determining significance:
  groupCoefs$z <- abs(groupCoefs$mean / groupCoefs$se) # z-value to evaluate significance
  ## Compute textual label:
  groupCoefs$zLabel <- as.character(round(groupCoefs$mean, nRound)) # copy over, to string
  ## Pad to desired string length:
  groupCoefs$zLabel <- ifelse(groupCoefs$mean > 0, str_pad(groupCoefs$zLabel, width = nRound + 2, pad = "0", side = "right"), groupCoefs$zLabel) # pad to nRound digits after comma
  groupCoefs$zLabel <- ifelse(groupCoefs$mean < 0, str_pad(groupCoefs$zLabel, width = nRound + 3, pad = "0", side = "right"), groupCoefs$zLabel) # pad to nRound digits after comma
  ## Handle cases of zero (no dot):
  groupCoefs$zLabel <- ifelse(!grepl("\\.", groupCoefs$zLabel), paste0("0.", paste0(rep("0", nRound), collapse = "")), groupCoefs$zLabel)
  ## Add stars and crosses for significance: 
  groupCoefs$zLabel <- ifelse(groupCoefs$z > 1.64 & groupCoefs$z < 1.96, paste0(groupCoefs$zLabel, "\U207A"), groupCoefs$zLabel) # latin cross
  groupCoefs$zLabel <- ifelse(groupCoefs$z > 1.96, paste0(groupCoefs$zLabel, "*"), groupCoefs$zLabel)
  groupCoefs$zLabel <- ifelse(groupCoefs$z > 3, paste0(groupCoefs$zLabel, "*"), groupCoefs$zLabel)
  
  # Alternatives to cross for marginally significant effects:
  # https://unicode-table.com/en/sets/crosses/
  # for (iLabel in 1:nrow(groupCoefs)){
  #   if (groupCoefs$z[iLabel] > 1.64 & groupCoefs$z[iLabel] < 1.96){
  #     groupCoefs$zLabel[iLabel] <- expression(paste(eval(groupCoefs$zLabel[iLabel]), "^+"))
  #     # groupCoefs$zLabel[iLabel] <- substitute(expression(n "^+"), list(n = groupCoefs$zLabel[iLabel]))
  #     # groupCoefs$zLabel[iLabel] <- bquote(.(groupCoefs$zLabel[iLabel])^+)
  #   }
  # } 
  
  # substitute(expression(a + b), list(a = groupCoefs$zLabel[iLabel]))
  # substitute(expression(a ^+), list(a = groupCoefs$zLabel[iLabel]))
  
  # groupCoefs$zLabel <- ifelse(groupCoefs$z > 1.64 & groupCoefs$z < 1.96, expression(paste0(groupCoefs$zLabel, "^+")), groupCoefs$zLabel) # latin cross
  # groupCoefs$zLabel <- ifelse(groupCoefs$z > 1.64 & groupCoefs$z < 1.96, paste0(groupCoefs$zLabel, "\U2670"), groupCoefs$zLabel) # latin cross
  # groupCoefs$zLabel <- ifelse(groupCoefs$z > 1.64 & groupCoefs$z < 1.96, paste0(groupCoefs$zLabel, "\U207A"), groupCoefs$zLabel) # superscript + (rather small)
  
  ## Drop intercept or not:
  if(dropIntercept){
    groupCoefs <- groupCoefs[2:nrow(groupCoefs), ] 
    # labels <- labels[2:length(labels)]
  }
  
  ## Determine final number of effects:
  nEff <- nrow(groupCoefs)
  
  ## Adjust number of colors:
  if (length(selCol) > nEff){selCol <-  selCol[1:nEff]} # if too many colors: only use first
  selCol <- rep(selCol, length.out = nEff) # if too few colors: repeat until enough
  # selCol <- rev(COL2("RdBu", nEff))
  
  ## Reverse order (first regressor will be plotted on top):
  if (revOrder){
    groupCoefs <- groupCoefs[nrow(groupCoefs):1,]
    selCol <- rev(selCol)
  }
  
  ## Compute index and lower/upper confidence bound:
  groupCoefs$idx <- seq(1, nrow(groupCoefs), 1) # numerical index of each effect to loop through (in correct order)
  groupCoefs$lower <- groupCoefs$mean - groupCoefs$se * SEweight
  groupCoefs$upper <- groupCoefs$mean + groupCoefs$se * SEweight
  
  # -------------------------------------------------------------------------- #
  ## Group-level plot:
  
  p <- ggplot(groupCoefs, aes(x = mean, y = label)) # define ggplot and axed
  
  ## A) Add error bar lines:
  cat("Plot error bar whiskers\n")
  for (iEff in 1:nEff){ # iEff <- 1
    ## For this effect: extract index, upper and lower end of whisker
    effData <- data.frame(x = c(groupCoefs$lower[iEff], groupCoefs$upper[iEff]),
                          y = rep(groupCoefs$idx[iEff], 2))
    p <- p + geom_line(data = effData, aes(x = x, y = y), size = 1.2, color = selCol[iEff])
  }
  
  # -------------------------------------------------------------------------- #
  ## B) Add fixed-effects points:
  
  cat("Plot fixed-effect coefficients\n")
  p <- p + geom_point(aes(x = mean, y = idx, color = factor(idx)), size = groupDotSize) +  # points for point estimates; size = 5
    scale_color_manual(values = selCol)
  
  # -------------------------------------------------------------------------- #
  ## C) Add subject-level points:
  
  ## Drop intercept or not (do outside plotSub for xLim determination):
  if(dropIntercept){
    subCoefs <- subCoefs[, 2:ncol(subCoefs)] 
    # labels <- labels[2:length(labels)]
  }
  if (is.vector(subCoefs)){subCoefs <- data.frame(subCoefs)} # ensure it has a column dimension
  
  ## Reverse order (first regressor will be plotted on top):
  if (revOrder){
    subCoefs <- subCoefs[, ncol(subCoefs):1]
  }
  if (is.vector(subCoefs)){subCoefs <- data.frame(subCoefs)} # ensure it has a column dimension
  
  if (plotSub){
    
    cat("Plot effect per subject\n")
    for (iEff in 1:nEff){ # iEff <- 1
      ## For this effect: extract and plot effects per subject
      effData <- data.frame(x = subCoefs[, iEff], # per-subject effect
                            y = rep(groupCoefs$idx[iEff], nSub) + subDisplacement) # y-axis position
      effData$y <- jitter(effData$y, amount = subJitter) # add jitter to distinuish subjects
      p <- p + geom_point(data = effData, aes(x = x, y = y), size = 2, 
                          # shape = 16, color = "gray30", # all grey
                          # shape = 21, color = "black", fill = "gray70", stroke = 1.2, # black edge, white fill
                          shape = 1, color = "black", # or color = selCol[iEff],
                          alpha = colAlpha)
    }  
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Determine x-axis dimensions for scaling:
  
  ## Extract range of subject-level effects for axis limits:
  if (!is.null(xLim)){ # if provided: retrieve from input
    xMin <- xLim[1]
    xMax <- xLim[2]
  } else { # else determine empirically based on subject coefficients
    ## Extract:
    if (plotSub){
      xMin <- min(subCoefs, na.rm = T)
      xMax <- max(subCoefs, na.rm = T)
    } else {
      xMin <- min(groupCoefs$lower, na.rm = T)
      xMax <- max(groupCoefs$upper, na.rm = T)
    }
    ## Erode a bit (need 10% for printing text for most positive coefficient):
    if(xMin > 0){xMin <- xMin * 0.90} else (xMin <- xMin * 1.10)
    if(xMax > 0){xMax <- xMax * 1.10} else (xMax <- xMax * 0.90)
    ## Assign:
    xLim <- c(xMin, xMax)
  }
  
  # -------------------------------------------------------------------------- #
  ## D) Add group-level coefficient as text:
  
  if (plotText){
    textDisplacement <- (xMax - xMin) * 0.10 # 10% of x-axis width
    cat("Print values of group-level effects as text\n")
    p <- p + geom_text(data = groupCoefs, aes(x = mean, y = idx, label = zLabel),  
                       nudge_x = textDisplacement, nudge_y = textOffset, na.rm = T, check_overlap = T, size = FTS/3) # nudge_x = 0.20
  }
  
  # -------------------------------------------------------------------------- #
  ### Other settings:
  
  ## Horizontal line at x = 0:  
  p <- p + geom_vline(xintercept = 0, 
                      linetype = "dashed", colour = "#949494", lwd = 1) # line at zero
  
  # -------------------------------------------------------------------------- #
  ## X-axis ticks:
  
  # xMin = -0.049; xMax = 0.459
  cat(paste0("xMin = ", round(xMin, 3), "; xMax = ", round(xMax, 3), "\n"))
  xRange <- xMax - xMin # determine range
  
  ## If symmetric: 5 symmetric break points
  if (abs(xMin) == abs(xMax)){
    
    breakVec <- sort(c(xMin, (xMin + mean(xLim))/2, mean(xLim), (mean(xLim) + xMax)/2, xMax)) # extremes and middle between extreme and center
    breakVec <- round(breakVec, 2) # round to 2 decimals
    
    ## else: determine adaptively:
  } else {
    
    ## Determine x-axis lower limit:
    iMag <- log10(abs(xMin)) # determine order of magnitude for rounding
    iMag <- ifelse(iMag < 0, floor(iMag), ceil(iMag)) # round
    xUnit <- 10 ^ iMag
    xBreakMin <- floor(xMin/xUnit)*xUnit # remove post-digit part, round, transform back
    
    ## Determine x-axis upper limit:
    iMag <- log10(abs(xMax)) # determine order of magnitude for rounding
    iMag <- ifelse(iMag < 0, floor(iMag), ceil(iMag)) # round
    xUnit <- 10 ^ iMag
    xBreakMax <- ceiling(xMax/xUnit)*xUnit # remove post-digit part, round, transform back
    
    ## xStep either 1 or 5; try out both:
    expVec <- seq(-3, 3, 1) # exponents for candidate ticks
    xTickVec <- sort(c(10^expVec, 5 * 10^expVec)) # candidate xSteps: either 1 or 5 of different magnitudes
    nTickTarget <- 4 # number desired ticks on x-axis
    xStepTarget <- (xBreakMax - xBreakMin) / nTickTarget # desired length xStep
    xStepIdx <- which(abs(xTickVec - xStepTarget) == min(abs(xTickVec - xStepTarget))) # find minimum
    xStepIdx <- xStepIdx[1] # first in case of multiple minima
    xStep <- xTickVec[xStepIdx] # extract optimal xStep
    
    ## Correct if one limit smaller than xStep:
    if (abs(xBreakMin) < xStep){xBreakMin <- 0}
    if (abs(xBreakMax) < xStep){xBreakMax <- 0}
    
    cat(paste0("xBreakMin = ", xBreakMin, ", xBreakMax = ", xBreakMax, ", xStep = ", xStep, "\n"))
    
    ## Distance between x-axis ticks:
    breakVec <- seq(xBreakMin, xBreakMax, xStep) # just very broad break points, aligned to magnitude
    
  }
  
  ## Axes:
  p <- p + coord_cartesian(xlim = xLim) # set x limits
  p <- p + scale_x_continuous(breaks = breakVec)
  p <- p + scale_y_continuous(breaks = 1:nEff, labels = groupCoefs$label)
  
  ## Labels:
  p <- p + labs(x = xLab,
                y = yLab)
  
  ## Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  ## Theme:
  p <- p + theme_classic() # base_size = 14
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(colour = "black", size = FTS),
                 axis.line = element_line(colour = "black"), # , linewidth = LWD), # fixed font sizes
                 axis.title = element_text(colour = "black", size = FTS), 
                 plot.title = element_text(colour = "black", size = FTS, hjust = 0.5), # center title 
                 legend.position = "none")
  
  print(p)
  cat("Finished :-)\n")
  
  return(p)
  
}

# ============================================================================ #
#### Plot intercorrelations of regressors in design matrix: #####

corrplot_regressors <- function(mod, perSub = F, varNames = NULL, savePNG = TRUE){
  #' Plot intercorrelations between regressors in design matrix using coefplot.
  #' @param mod model object fitted with lme4 or another package.
  #' @param perSub compute correlation between regressors separately per subject, than average (TRUE) (default: FALSE).
  #' @param varNames vector of strings, names of variables to use for rows/ columns of design matrix.
  #' @param savePNG Boolean, save as PNG to dirs$plot (TRUE, default) or not (FALSE).
  #' @return nothing, plots to console and saves in dirs$plot.
  
  require(psych) # for fisherz and fisherz2r
  require(corrplot)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Extract design matrix:
  
  DM <- model.matrix(mod) # extract design matrix
  DM <- DM[, 2:ncol(DM)] # drop intercept
  
  # -------------------------------------------------------------------------- #
  ## Compute correlation:
  
  if (perSub){
    
    cat("Separate DM per subject, average\n")
    subIdx <- mod@frame[, ncol(mod@frame)] # extract subject indices
    stopifnot(nrow(DM) == length(subIdx)) # assure dimensions match
    subVec <- unique(subIdx)
    nSub <- length(subVec)
    
    Mlist <- vector(mode = "list", length = nSub) # initialize
    
    for (iSub in 1:nSub){ # iSub <- 1
      
      subID <- subVec[iSub] # subject ID
      DMsub <- DM[subIdx == subID, ] # select part of design matrix for this subject
      M <- cor(DMsub) # correlation
      MF <- fisherz(M) # Fisher-z transform
      diagIdx <- which(MF == Inf) # detect diagonal elements (infinite)
      MF[diagIdx] <- 1 # temporarily overwrite to 1, correct later after transforming back
      
      Mlist[[iSub]] <- MF # store
    }
    
    M <- Reduce("+", Mlist)/nSub # mean across subjects
    M <- fisherz2r(M) # transform back
    M[diagIdx] <- 1 # set diagonal back to 1 
    
  } else {
    M <- cor(DM)
  }
  
  ## Print range to console:
  Mvec <- as.vector(M) # to vector
  diagVec <- seq(1, length(Mvec), nrow(M) + 1) # identify indices of diagonal
  Mvec[diagVec] <- NA # set diagonal to NA
  nRound <- 2
  cat(paste0("All correlations between r = ", round(min(Mvec, na.rm = T), nRound), " and r = ", round(max(Mvec, na.rm = T), nRound), "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Overwrite variables names:
  
  if (is.null(varNames)){
    rownames(M) <- substitute_label(rownames(M)) # substitute for known variable names
  } else {
    stopifnot(nrow(M) == length(varNames)) # check if same length
    rowNames(M) <- varNames # overwrite
  }
  colnames(M) <- rownames(M)
  
  # -------------------------------------------------------------------------- #
  ## Title and name for saving:
  
  # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
  if (class(mod) %in% c("glmerMod")){
    
    formulaStr <- mod@call$formula
    
  } else {
    
    formulaStr <- attr(mod@frame, "formula")
    
  }
  formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
  
  # https://stackoverflow.com/questions/40509217/how-to-have-r-corrplot-title-position-correct
  # titleStr <- paste0("Intercorrelation regressors for \n", deparse(formulaStr), width.cutoff = 20))
  
  plotName <- paste0("corr_reg_", formula2handle(formulaStr))
  if (perSub){plotName <- paste0(plotName, "_perSub")}
  plotName <- paste0(plotName, ".png")
  
  # -------------------------------------------------------------------------- #
  ## Make corrplot:
  
  # mar = c(0,0, 1,0)  
  # if(savePNG) {png(paste0(dirs$plot, "ic_reg/", plotName), width = 480, height = 480)}
  if(savePNG) {png(paste0(dirs$plot, plotName), width = 480, height = 480)}
  
  # https://stackoverflow.com/questions/40352503/change-text-color-in-corrplot-mixed
  
  # corrplot(M, method = "circle", col = rev(COL2('RdBu', 200))) # colored dots of different size
  # corrplot(M, method = "number", col = rev(COL2('RdBu', 200))) # numerals of different color
  
  ## Upper half colored dots of different size, lower half black numerals, variable names in diagonal:
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = "black", upper.col = rev(COL2('RdBu', 200)), tl.col = "black")
  
  ## Colors dots of different size with black numerals in them:
  corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  # corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.offset = 1, tl.srt = 0) # column labels higher, not rotated
  
  ## Also numerals in color, different color scale (uniform), variable names in diagonal:
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = COL1('YlOrRd', 200), upper.col = rev(COL2('RdBu')))
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = COL1('YlOrRd', 200), upper.col = COL1('YlOrRd', 200))
  
  # print(p)
  if(savePNG){
    dev.off(); 
    cat(paste0("Saved under \n", plotName, " :-)\n"))
    corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  }
  
}

# ============================================================================ #
#### Plot intercorrelations of coefficients from model: #####

corrplot_coefficients <- function(input, varNames = NULL, savePNG = TRUE){ 
  #' Plot intercorrelations between regressors in design matrix using coefplot.
  #' @param mod model object fitted with lme4 or another package.
  #' @param varNames vector of strings, names of variables to use for rows/ columns of design matrix.
  #' @param savePNG Boolean, save as PNG to dirs$plot (TRUE, default) or not (FALSE).
  #' @return nothing, plots to console and saves in dirs$plot.
  
  require(corrplot)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ### Detect class and extract coefficients:
  
  ## Detect model class:  
  modClass <- class(input)
  if (modClass == "list"){modClass <- class(input$modList[[1]]); modClass <- modClass[1]}
  cat(paste0("Input model of class ", modClass, "\n"))
  
  ## Extract coefficients, parameter names, formula:
  
  if(modClass %in% c("glmerMod", "lmerMod", "lmerTest", "lmerModLmerTest")){ # from lme4
    
    coefMat <- coef(input)[[1]]
    
    parNamesVec <- colnames(coefMat)
    
    if (modClass %in% c("glmerMod")){
      
      formulaStr <- input@call$formula
      
    } else {
      
      formulaStr <- attr(input@frame, "formula")
      
    }
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    # lmer
  } else if (modClass %in% c("mixed")){ # from afex
    
    coefMat <- coef(input$full_model)[[1]]
    parNamesVec <- rownames(coefMat)
    
    formulaStr <- input@call$formula
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
  } else if (modClass %in% c("lm", "glm")){ # from lm
    
    coefMat <- input$bMat
    parNamesVec <- names(coef(input$modList[[1]]))
    if (modClass == "glm"){
      formulaStr <- input$modList[[1]]$formula
    } else if (modClass == "lm"){
      formulaStr <- eval(input$modList[[1]]$call[[2]]) 
    } else {
      stop("Unknown model class")
    }
    
  } else if (modClass %in% c("brmsfit")){ # from brms
    
    ## Parameter names and formula:
    parNamesVec <- row.names(fixef(input)) # names of all predictors
    nParam <- length(parNamesVec)
    
    formulaStr <- input$formula
    formulaStr <- formulaStr[[1]] # extract only first object
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    ## Exttract correlation estimates:
    brmsVarCorr <- VarCorr(input)[[1]]$cor 
    brmsVarCorr <- VarCorr(input)$subject_f$cor 
    
    coefMat <- matrix(NA, nParam, nParam) # initialize
    rownames(coefMat) <- parNamesVec
    colnames(coefMat) <- parNamesVec
    for (iParam1 in 1:nParam){ # iParam1 <- 1
      for (iParam2 in 1:nParam){ # iParam2 <- 2
        coefMat[iParam1, iParam2] <- brmsVarCorr[iParam1, 1, iParam2]
      }
    }
    
  } else {
    
    stop("Unknown model class")
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Compute correlation:
  
  M <- cor(coefMat)
  
  ## Print range to console:
  Mvec <- as.vector(M) # to vector
  diagVec <- seq(1, length(Mvec), nrow(M) + 1) # identify indices of diagonal
  Mvec[diagVec] <- NA # set diagonal to NA
  nRound <- 2
  cat(paste0("All correlations between r = ", round(min(Mvec, na.rm = T), nRound), " and r = ", round(max(Mvec, na.rm = T), nRound), "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Overwrite variables names:
  
  if (is.null(varNames)){
    
    rownames(M) <- parNamesVec
    rownames(M) <- substitute_label(rownames(M)) # substitute for known variable names
    
  } else { # overwrite with inputs
    
    stopifnot(nrow(M) == length(varNames)) # check if same length
    rowNames(M) <- varNames # overwrite
  }
  
  colnames(M) <- rownames(M) # copy over
  
  # -------------------------------------------------------------------------- #
  ## Title and name for saving:
  
  # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
  
  plotName <- paste0("corr_coef_", modClass, "_", formula2handle(formulaStr), ".png")
  # plotNameFull <- paste0(dirs$plot, "ic_coef/", plotName)
  plotNameFull <- paste0(dirs$plot, plotName)
  cat(paste0("File path has ", nchar(plotNameFull), " characters\n"))
  if (nchar(plotNameFull) > 260){
    warning("File path too long, shorten\n")
    plotNameFull <- paste0(substr(plotNameFull, 1, 255)[1], ".png")
  }
  
  # -------------------------------------------------------------------------- #
  ## Visualize correlation matrix:
  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  
  if(savePNG) {png(plotNameFull, width = 480, height = 480)}
  
  # corrplot(M, method = "circle", col = rev(COL2('RdBu', 200)))
  # corrplot(M, method = "number", col = rev(COL2('RdBu', 200)))
  corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  
  if(savePNG){
    dev.off(); 
    cat(paste0("Saved under \n", plotName, " :-)\n"))
    ## Plot again:
    corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  }
  
}

# ============================================================================ #
#### Save coefficients from model: #####

save_coefficients <- function(input){ 
  #' Save coefficients from model output as .csv-file.
  #' @param mod model object fitted with lme4 or another package.
  #' @return nothing, plots to console and saves in dirs$plot.
  
  # -------------------------------------------------------------------------- #
  ### Detect class and extract coefficients:
  
  ## Detect model class:  
  modClass <- class(input)
  if (modClass == "list"){modClass <- class(input$modList[[1]]); modClass <- modClass[1]}
  cat(paste0("Input model of class ", modClass, "\n"))
  
  ## Extract coefficients and formula:
  
  if(modClass %in% c("glmerMod", "lmerMod", "lmerTest", "lmerModLmerTest")){ # from lme4
    
    coefMat <- coef(input)[[1]]
    
    if (modClass %in% c("glmerMod")){
      
      formulaStr <- input@call$formula
      
    } else {
      
      formulaStr <- attr(input@frame, "formula")
      
    }
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    # lmer
  } else if (modClass %in% c("mixed")){ # from afex
    
    coefMat <- coef(input$full_model)[[1]]
    
    formulaStr <- input@call$formula
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
  } else if (modClass %in% c("lm", "glm")){ # from lm
    
    coefMat <- input$bMat
    colnames(coefMat) <- names(coef(input$modList[[1]])) # add variable names
    if (modClass == "glm"){
      formulaStr <- input$modList[[1]]$formula
    } else if (modClass == "lm"){
      formulaStr <- eval(input$modList[[1]]$call[[2]]) 
    } else {
      stop("Unknown model class")
    }
    
  } else if (modClass %in% c("brmsfit")){ # from brms
    
    ## Parameter names and formula:
    
    formulaStr <- input$formula
    formulaStr <- formulaStr[[1]] # extract only first object
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    ## Exttract mean of per-subject posterior:
    coefArray <- coef(input)[[1]] # extract array of dimensions nSub x 4 x nParam
    nSub <- dim(coefArray)[1]
    nParam <- dim(coefArray)[3]
    
    coefMat <- matrix(NA, nSub, nParam) # initialize
    colnames(coefMat) <- row.names(fixef(input))
    for (iSub in 1:nSub){ # iParam1 <- 1
      for (iParam in 1:nParam){ # iParam2 <- 2
        coefMat[iSub, iParam] <- coefArray[iSub, 1, iParam]
      }
    }
    
  } else {
    
    stop("Unknown model class")
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Name for saving:
  
  fileName <- paste0("coefMat_", modClass, "_", formula2handle(formulaStr), ".csv")
  cat(paste0("Save data under ", fileName, "\n"))
  fileNameFull <- paste0(dirs$outputDir, fileName)
  cat(paste0("File path has ", nchar(fileNameFull), " characters\n"))
  if (nchar(fileNameFull) > 260){
    warning("File path too long, shorten\n")
    fileNameFull <- paste0(substr(fileNameFull, 1, 255)[1], ".csv")
  }
  
  # -------------------------------------------------------------------------- #
  ## Save:
  
  write.csv(coefMat, fileNameFull, row.names = F)
  
  cat("Finished :-)\n")
  
}

# ============================================================================ #
#### Find plausible tick step size given axis limits: #####

find_step <- function(xLim, nTickTarget = 5){
  #' Plot learning curve per cue per subject using basic plot() function 
  #' @param xLim vector of 2 numerics, axis limits.
  #' @param nTickTarget scalar integer, number of desired axis ticks (default: 5).
  #' @return xStep scalar numeric, optimal step size
  
  if(length(xLim) != 2){stop("xLim must have 2 elements")}
  if(length(xLim) != 2){stop("nTickTarget must have 1 element")}
  if(nTickTarget != round(nTickTarget)){stop("nTickTarget must be an integer")}
  
  expVec <- seq(-5, 5, 1) # exponents for candidate ticks
  xTickVec <- sort(c(10^expVec, 5 * 10^expVec)) # candidate xSteps: either 1 or 5 of different magnitudes
  xStepTarget <- (xLim[2] - xLim[1]) / nTickTarget # desired length xStep
  xStepIdx <- which(abs(xTickVec - xStepTarget) == min(abs(xTickVec - xStepTarget))) # find minimum
  xStepIdx <- xStepIdx[1] # first in case of multiple minima
  xStep <- xTickVec[xStepIdx] # extract optimal xStep
  return(xStep)
  
}

# ============================================================================ #
#### Find plausible axis limits: #####

find_round_lim <- function(input){
  #' Round input xMin and xMax such that they become plausible axis limits 
  #' given their magnitude. 
  #' @param input vector of 2 numerics, axis limits.
  #' @return xStep scalar numeric, optimal step size
  
  # input
  signVec <- ifelse(input > 0, 1, -1)
  absVec <- abs(input)
  xUnit <- min(floor(log10(absVec))) # detect smaller of both magnitudes
  xScale <- 10^xUnit # scaling factor 
  output <- round(input/xScale)*xScale # scale up/down, round, scale down/up again
  
  while(output[1] == output[2]){
    xUnit <- xUnit - 1 # reduce unit 
    xScale <- 10^xUnit # scaling factor 
    output <- round(input/xScale)*xScale # scale up/down, round, scale down/up again
  }
  
  return(output)
}

# =============================================================================================== #
#### Quick CIs based on SEs: ####

quickCI <- function(mod, selEff = NULL, level = 0.95, nRound = 2){
  #' Compute CIs for given lme4 model given SEs from model.
  #' @param data mod model objected fitted with lme4.
  #' @param selEff vector of integers, index of effect in model for which to compute effect size (default: 2).
  #' @param savePNG numeric, 0-1, level of CIs (default: 0.95).
  #' @param nRound integer, number of digits after comma to round to (default: 2).
  #' @return print to console.
  require(arm) # for se.fixef
  
  twoSideLevel <- 1 - (1 - level) / 2 # correct for two-sided test
  zVal <- qnorm(twoSideLevel) # respective z-level threshold
  
  ## If no effect selected: print all effects in model (skipping intercept)
  if (is.null(selEff)){selEff <- 2:length(fixef(mod))}
  
  ## Loop through effects:
  for (iEff in selEff){
    
    # print(round(c(fixef(mod)[iEff] - zVal*se.fixef(mod)[iEff], fixef(mod)[iEff] + zVal*se.fixef(mod)[iEff]), nRound))
    tmp <- round(c(fixef(mod)[iEff] - zVal*se.fixef(mod)[iEff], 
                   fixef(mod)[iEff] + zVal*se.fixef(mod)[iEff]), 
                 nRound)  
    cat(paste0(level*100, "%-CIs for ", colnames(model.matrix(mod))[iEff], 
               ": b = ", round(fixef(mod)[iEff], 3), ", ", 
               level*100, "%-CI [", paste(tmp, collapse = ", "), "]\n"))
  }
}

# =============================================================================================== #
#### Print effect from lme4 model: #####

print_effect <- function(mod, eff, nDigit = 3){
  #' Print selected effect from lme4 model
  #' @param mod fitted model
  #' @param eff string, name of effect for which to print effect
  #' @param nDigit integer, number of digits to round after comma, default 2
  #' @return nothing returned, but printed
  require(stringr)
  
  nPad <- nDigit + 2
  
  if (str_sub(eff,-1) == "f"){eff <- paste0(eff, "1")} # add 1 at the end if eff is factor
  
  # Extract output of fixed effects:
  coefs <- summary(mod)$coefficients # extract coefficients
  idx <- which(rownames(coefs)==eff) # find effect back
  if (length(idx)==0){stop(paste0("Effect ", eff, " not found"))} 
  
  ## Retrieve coefficients:
  if (summary(mod)$objClass== "glmerMod"){ # glmer
    
    # Extract relevant info:
    b <- coefs[idx, 1]
    se <- coefs[idx, 2]
    zScore <- coefs[idx, 3]
    pVal <- coefs[idx, 4]
    
  } else if (summary(mod)$objClass== "lmerModLmerTest"){ # lmer
    
    # Extract relevant info:
    b <- coefs[idx, 1]
    se <- coefs[idx, 2]
    dfs <- coefs[idx, 3]
    zScore <- coefs[idx, 4]
    pVal <- coefs[idx, 5]
    
  } else {
    stop("Unknown model type")
  }
  
  # Variable padding of b based on sign:
  bPad <- ifelse(b > 0, nPad, nPad+1) # pad to 5 digits if negative
  zPad <- ifelse(zScore > 0, nPad, nPad+1) # pad to 5 digits if negative
  
  ## Handle b:
  if (round(b, nDigit) == 0){
    bText <- "0"
  } else {
    bText <- str_pad(round(b, nDigit), bPad, side = "right", pad = "0")
  }
  
  ## Handle se:
  if (round(se, nDigit) == 0){
    seText <- "0"
  } else {
    seText <- str_pad(round(se, nDigit), nPad, side = "right", pad = "0")
  }
  
  ## Handle statistic for given object:
  if (summary(mod)$objClass == "glmerMod"){
    zStat <- ", z = "
  } else {
    zStat <- paste0(", t(",round(dfs, nDigit), ") = ")
  }
  
  if (round(zScore, nDigit) == 0){
    zText <- "0"
  } else {
    zText <- str_pad(round(zScore, nDigit), zPad, side = "right", pad = "0")
  }
  
  ## Handle very small p-values:
  if (pVal < 0.001){
    pText <- "p < .001"
  } else {
    pText <- paste0("p = ", str_pad(round(pVal,(nDigit+1)), 5, side = "right", pad = "0")) # p-value: always 5 digits
  }
  
  # Print to console:
  cat(paste0("b = ", bText,
             ", se = ", seText,
             zStat, zText,
             ", ", pText, "\n"))
}

# ============================================================================ #
#### Fit lm per subject: #####

loop_lm_subject <- function(data, formula, isBinomial = F, family = "binomial",
                            subVar = "subject_n"){
  #' Perform lm separately for each subject, store coefficients and models, 
  #' one-sample t-test across subjects for each effect, return.
  #' @param data data frame with variable subject and DVs and IVs.
  #' @param formula string with formula to fit in Wilkinson notation.
  #' @param isBinomial boolean, fit generalized lm with binomial link function (T) or not (F).
  #' @param family distribution of DV to use (default: binomial).
  #' @return output list with elements:
  #' output$bVec: vector of size nSub x nEff, b weights for each effect for each subject.
  #' output$modList: list of models of each subject.
  #' Prints results from t-test across subjects for each effect.
  
  require(DescTools)
  
  # ----------------------------------------------- #
  ## Fixed variables:
  if (!(subVar %in% names(data))){stop(paste0("subVar ", subVar, "not contained in data"))}
  nDigit <- 3
  
  # ----------------------------------------------- #
  ## Determine number of subjects:
  subVec <- unique(data[, subVar])
  nSub <- length(subVec)
  cat(paste0("Found data from ", nSub, " subjects\n"))
  
  # ----------------------------------------------- #
  ## Fit model for first subject available to determine number of coefficients:
  subIdx <- which(data[, subVar] == subVec[1])
  subData <- data[subIdx, ]
  mod <- lm(formula = formula, data = subData) # fit model
  nEff <- length(mod$coefficients)
  
  # ----------------------------------------------- #
  ## Initialize matrix to store coefficients:
  
  bMat <- matrix(NA, nrow = nSub, ncol = nEff) # initialize
  modList <- list()
  sumModList <- list()
  
  ## Loop through all subjects:
  for (iSub in 1:nSub){ # iSub <- 1
    
    ## Select subject and data:    
    subID <- subVec[iSub]
    cat(paste0("Start subject ", subID, "\n"))
    subIdx <- which(data[, subVar] == subID)
    subData <- data[subIdx, ]
    
    ## Fit model:
    if (isBinomial) {
      if (family == "binomial"){
        mod <- glm(formula = formula, data = subData, family = binomial())
      }
      else if (family == "poisson"){
        mod <- glm(formula = formula, data = subData, family = poisson())
      } else {
        stop("Unknown family for DV")
      }
    } else {
      mod <- lm(formula = formula, data = subData)
    }
    
    modList[[iSub]] <- mod
    sumModList[[iSub]] <- summary(mod)
    bMat[iSub, ] <- mod$coefficients # store data
    
  } # end iSub
  
  # names(bMat) <- names(coef(mod)) # copy over regressor names
  
  # ----------------------------------------------- #
  ## Perform 1-sample t-test across subjects for each effect:
  
  if(nEff > 0) {
    for (iEff in 1:nEff){ # iEff <- 1
      # out <- t.test(FisherZ(bMat[, iEff]) # one-sample t-test: only for correlations in range [-1, 1], so not for intercept, not for glm
      out <- t.test((bMat[, iEff])) # one-sample t-test
      cat(paste0("Effect for ", names(mod$coefficients)[iEff], 
                 ": t(", out$parameter, ") = ", round(out$statistic, nDigit), 
                 ", p = ", out$p.value, "\n"))
    }
    
  } else {
    cat("Intercept-only model, no effects printed\n")
  }
  
  cat("Finished! :-)\n")
  
  # ----------------------------------------------- #
  ## Output:
  output <- list()
  output$bMat <- bMat
  output$modList <- modList
  output$sumModList <- sumModList
  
  return(output)
}

# =============================================================================================== #
#### Fit lm per subject: #####

loop_glm_subject <- function(data, formula, family = "binomial", subVar = "subject_n"){
  #' Wrapper on loop_lm_subject to perform glm separately for each subject.
  #' @param data data frame with variable subject and DVs and IVs.
  #' @param formula string with formula to fit in Wilkinson notation.
  #' @param family distribution of DV to use (default: binomial).
  #' @return modList list with elements:
  #' modList$bVec: vector of size nSub x nEff, b weights for each effect for each subject.
  #' modList$modList: list of models of each subject.
  #' Prints results from t-test across subjects for each effect.
  
  ## Call loop_lm_subject with isGLM = T:
  out <- loop_lm_subject(data, formula, isBinomial = T, family = family, subVar = subVar)
  
  return(out)
}

# =============================================================================================== #
#### Aggregate data per subject per condition: ####

aggregate_sub_cond <- function(data, yVar, xVarVec, subVar = "subject_f",
                               format = "long", printSumStats = F){
  #' Compute mean and SD of yVar per conditions spanned by xVarVec.
  #' @param data data frame, input trial-by-trial data.
  #' @param yVar scalar string, name of dependent variable which to aggregate.
  #' @param xVarVec vector of strings, names of independent variables by which to aggregate.
  #' @param subVar scalar string, name of subject identifier (default: subject_f).
  #' @return nothing, just print to console
  require(plyr)
  
  # yVar <- "RTcleaned_n"
  # xVarVec <- c("arousal_f", "reqAction_f", "valence_f")
  # subVar <- "subject_f"
  
  ## Copy over:
  input <- data
  
  ## Aggregate into long format:
  cat("Aggregate data per subject per condition into long format\n")
  cat(paste0("Conditions are formed by ", paste0(xVarVec, collapse = ", "), "\n"))
  aggrData_long <- eval(parse(text = paste0("ddply(input, .(subject_f, ", paste0(xVarVec, collapse = ", "), "), function(subData){
    ", yVar, " <- mean(subData[, \"", yVar, "\"], na.rm = T)
    return(data.frame(", yVar, "))
    dev.off()})")))
  output <- aggrData_long
  
  ## Create overall condition variable by concenating xVars:
  if (format == "wide" | format == "WIDE"){
    cat("Cast into wide format\n")
    aggrData_long$condition_f <- eval(parse(text = paste0("paste0(aggrData_long$", paste0(xVarVec, collapse = ", \"_\", aggrData_long$"), ")")))
    aggrData_long <- aggrData_long[, c(subVar, yVar, "condition_f")]
    
    ## Reshape into wide format:
    aggrData_wide <- reshape(aggrData_long, direction = "wide",
                             idvar = subVar, v.names = yVar, timevar = "condition_f")
    output <- aggrData_wide
    
    if(printSumStats){
      ## Print condition means and SDs:
      nRound <- 3
      condData <- aggrData_wide[, grepl(yVar, names(aggrData_wide))]
      cat("Mean per condition:\n")
      # colMeans(condData, na.rm = T)
      print(round(apply(aggrData_wide[, grepl(yVar, names(aggrData_wide))], 2, mean, na.rm = T), nRound))
      cat("SD per condition:\n")
      print(round(apply(aggrData_wide[, grepl(yVar, names(aggrData_wide))], 2, sd, na.rm = T), nRound))
      cat("Finished! :-)\n")
      
    }
    
  }
  # str(aggrData_wide)
  
  return(output)
}

# =============================================================================================== #
#### Aggregate and run RM-ANOVA with ezANOVA: ####

run_RMANOVA <- function(data, yVar, xVarVec, zVar = NULL, subVar = "subject_f"){
  #' Remove NAs, aggregate per subject per condition, run RM-ANOVA with ez.
  #' 
  require(ez)
  
  input <- data # copy over
  
  aggrData_long <- aggregate_sub_cond(input, yVar = yVar, xVarVec = xVarVec, subVar = subVar, format = "long")
  
  ## Run RM-ANOVA:
  cat("Run RM-ANOVA\n")
  out <- eval(parse(text = paste0("ezANOVA(data = aggrData_long, dv = .(", yVar, "), within = .(", 
                                  paste0(xVarVec, collapse = ", "), "), wid = .(", subVar, "), type = 3, detailed = F)")))
  # , detailed = TRUE, type = \"III\")")))
  ## Print output to console:
  print(out$ANOVA)
  return(out$ANOVA)
}

# ================================================================================================================================================ #
#### Barplot 1 IV: Aggregate per condition per subject, plot (1 IV on x-axis): ####

custom_barplot1 <- function(data, xVar = NULL, yVar = NULL, subVar = "subject_n", 
                            xLab = NULL, yLab = NULL, main = NULL, selCol = NULL,
                            isPoint = T, isConnect = T, isMidLine = F, hLine = NULL, isBeeswarm = F, 
                            yLim = NULL, FTS = NULL, 
                            savePNG = T, saveEPS = F, prefix = NULL, suffix = NULL){
  #' Make bar plot with error bars and individual-subject data points.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: retrieve appropriate name with substitute_label()).
  #' @param yLab string, label for y-axis (default: retrieve appropriate name with substitute_label()).
  #' @param main string, overall plot label (optional).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: FALSE).
  #' @param isMidLine Boolean, add horizontal line at midpoint of y-axis (default: FALSE).
  #' @param yLine scalar numeric, draw horizontal line at given y value (default: NULL).   
  #' @param isConnect Boolean, connect individual data points with grey lines (default: FALSE).
  #' @param isBeewswarm Boolean, plot individual data points per condition as beeswarm densities (default: FALSE).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param FTS scalar integer, font size to use.
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @param prefix string, string to add at the beginning of plot name (optional).
  #' @param suffix string, string to add at the end of plot name (optional).
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load required packages:
  
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  ## Input variables:
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  
  ## Axis labels:
  if(is.null(xLab)){xLab <- substitute_label(xVar)}
  if(is.null(yLab)){yLab <- substitute_label(yVar)}
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  SEweight <- 1
  LWD <- retrieve_plot_defaults("LWD")
  dodgeVal <- retrieve_plot_defaults("dodgeVal")
  colAlpha <- 0.5 # 1
  
  if (is.null(FTS)){
    FTS <- retrieve_plot_defaults("FTS") # 30 or 15?
  }
  
  if (length(unique(data[, xVar])) > 20){
    cat("More than 20 levels on x-axis, reduce font size and line width\n")
    FTS <- 12
    LWD <- 0.5
  }
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  
  cat("Overall condition means (without first aggregating per subject):\n")
  print(tapply(data[, yVar], data[, xVar], mean, na.rm = T))
  
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$subject <- data[, subVar]
  
  ## Exclude NAs:
  nRow1 <- nrow(data)
  data <- droplevels(subset(data, !(is.na(x)) & !(is.na(y))))
  nRow2 <- nrow(data)
  cat(paste0("Excluded ", nRow1 - nRow2, " rows due to NAs\n"))
  
  ## Colours:
  if(is.null(selCol)){
    selCol <- retrieve_colour(xVar)
    selCol <- rep(selCol, length.out = length(unique(data[, xVar])))
  }
  if(length(selCol) != length(unique(data[, xVar]))){
    stop(paste0("Length selCol = ", length(selCol), " while number levels xVar = ", length(unique(data[, xVar])), ", do not match"))
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  
  aggrData <- ddply(data, .(subject, x), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  cat(paste0("Min = ", round(min(aggrData$y), 3), "; Max = ", round(max(aggrData$y), 3)), "\n")
  
  # -------------------------------------------------------------------------- #
  ## Add jittered x-axis variable for points:
  
  aggrData$xpos <- as.numeric(aggrData$x) # to numeric
  aggrData$xpos <- aggrData$xpos - min(aggrData$xpos) + 1 # plot starts at lowest level of x-variable
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = .09) # jitter 0.09
  aggrData$xj <- aggrData$xpos + aggrData$j # add jitter
  
  # -------------------------------------------------------------------------- #
  ## Determine y limits if not given:
  
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x"))
  
  # -------------------------------------------------------------------------- #
  ## Control settings for saving:
  
  ## Additions:
  if(is.null(prefix)){prefix <- ""} else {prefix <- paste0(prefix, "_")}
  if(is.null(suffix)){suffix <- ""} else {suffix <- paste0("_", suffix)}
  
  ## Name:
  plotName <- paste0("custombarplot1_", prefix, yVar, "_", xVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  plotName <- paste0(plotName, suffix)
  cat(paste0("Start plot ", plotName, "\n"))
  
  # Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  # -------------------------------------------------------------------------- #
  ## Start ggplot:
  p <- ggplot(summary_d,aes(x, y))
  
  ## Bars of means:
  cat("Add group-level bars \n")
  p <- p + stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.6,
                        lwd = LWD, fill = selCol, color = "black") + 
    
    ## Error bars:
    cat("Add error bars \n")
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se * SEweight, ymax = y + se * SEweight),
                         position = position_dodge(width = dodgeVal), width = 0.1,
                         lwd = LWD, color = "black", alpha = 1)
  
  # -------------------------------------------------------------------------- #
  ## Individual data points:
  
  cat("Add per-subject data points\n")
  if (isPoint){
    ## Colored dots:
    p <- p + geom_point(data = aggrData, aes(x = xj, fill = x), shape = 21, size = 2, stroke = 1.2, # size = 0.6, 
                        color = "black", alpha = colAlpha)
    p <- p + scale_fill_manual(values = selCol, limits = levels(aggrData$x))
    ## Grey dots:
    # p <- p + geom_point(data = aggrData, aes(x = xj), shape = 21, size = 2, fill = NA, stroke = 1.5, # size = 0.6, 
    #                     color = "grey", alpha = colAlpha) # color = black, grey60,
  }
  
  if (isConnect){
    ## Connect colored dots:
    
    subVec <- sort(unique(aggrData$subject))
    nSub <- length(subVec)
    for(iSub in 1:nSub){
      subData <- subset(aggrData, subject == subVec[iSub])
      p <- p + geom_path(data = subData, aes(x = xj, y= y), color = 'grey40', # color = 'grey70'
                         alpha = 0.40, lwd = 0.5) # alpha = 0.80, lwd = 1
    }
  }
  
  if (isBeeswarm){
    p <- p + geom_beeswarm(data = aggrData, aes(x = xpos), shape = 1, size = 2, stroke = 1, # size = 0.6, 
                           color = "black", alpha = colAlpha)
  }
  
  # -------------------------------------------------------------------------- #
  ### Settings:
  
  ## Add horizontal lines:
  if (isMidLine){
    yMid <- (yLim[1] + yLim[2])/2
    p <- p + geom_hline(yintercept = yMid, linetype = 2, color = "black", linewidth = 1) # Middle line at 0
  }
  
  if (!(is.null(hLine))){
    p <- p + geom_hline(yintercept = hLine, linetype = 2, color = "black")
  }
  
  ## Y-axis labels:
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  p <- p + coord_cartesian(ylim = yLim)
  
  ## Axis labels:
  p <- p + labs(x = xLab, y = yLab)
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  ## Theme:
  p <- p + theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 # axis.text.x = element_blank(), # remove x-axis labels
                 # axis.ticks.x = element_blank(), # remove x-axis labels
                 axis.title = element_text(size = FTS), 
                 title = element_text(size = FTS),
                 legend.position = "none",
                 axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
  print(p)
  if(savePNG | saveEPS){
    dev.off(); 
    print(p)
  }
  return(p)
  cat("Finished :-)\n")
}

# ================================================================================================================================================ #
#### Barplot 2 IVs: Aggregate per condition per subject, plot (2 IVs, 1 on x-axis, 1 as color): ####

custom_barplot2 <- function(data, xVar, yVar, zVar, subVar = "subject_n", 
                            xLab = NULL, yLab =  NULL, zLab = NULL, main = NULL,
                            selCol = NULL, 
                            isPoint = T, isConnect = T, isBeeswarm = F, addLegend = F,
                            yLim = NULL, FTS = NULL, dotSize = NULL,
                            savePNG = T, saveEPS = F, prefix = NULL, suffix = NULL){
  #' Make bar plots with 2 IVs: x-axis and color.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: retrieve appropriate name with substitute_label()).
  #' @param yLab string, label for y-axis (default: retrieve appropriate name with substitute_label()).
  #' @param zLab string, label for color legend (default: retrieve appropriate name with substitute_label()).
  #' @param main string, overall plot label (optional).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param isConnect Boolean, connect individual data points with grey lines (default: FALSE).
  #' @param isBeeswarm Boolean, plot individual data points per condition as beeswarm density (default: FALSE).
  #' @param addLegend Boolean, add legend for z-axis (colour) (default: TRUE).
  #' @param FTS scalar integer, font size to use (default: 30).
  #' @param dotSize scalar integer, size of single-subject dots to use (default: 1).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @param prefix string, string to add at the beginning of plot name (optional).
  #' @param suffix string, string to add at the end of plot name (optional).
  #' @return creates (and saves) plot.
  
  # xLab = NULL; yLab = NULL; zLab = NULL; main = NULL; selCol = NULL;
  # isPoint = T; isConnect = T; isBeeswarm = F; addLegend = TRUE; yLim = NULL; FTS = NULL; dotSize = NULL; savePNG = T; saveEPS = F; prefix = NULL; suffix = NULL
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  require(ggbeeswarm) # for ggbeeswarm
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  ## Input variables:
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(zVar %in% names(data))){stop("zVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  
  if(is.numeric(data[, xVar])){data[, xVar] <- factor(data[, xVar]); cat("Convert xVar to factor\n")}
  
  ## Axis labels:
  if(is.null(xLab)){xLab <- substitute_label(xVar)}
  if(is.null(yLab)){yLab <- substitute_label(yVar)}
  if(is.null(zLab)){zLab <- substitute_label(zVar)}
  zLab <- gsub(" ", " \n", zLab) # add newline to any z-label
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  LWD <- retrieve_plot_defaults("LWD") # 1.3
  if (is.null(FTS)){
    FTS <- retrieve_plot_defaults("FTS") # 30
    if (length(unique(data[, xVar])) > 10){FTS <- FTS/2; cat("Half font size because many x levels\n")}
  }
  if (is.null(dotSize)){
    dotSize <- retrieve_plot_defaults("dotSize") # 0.5
    if (length(unique(data[, xVar])) * length(unique(data[, zVar])) > 15){dotSize <- 0.5; cat("Lower dotSize from 1 to 0.5 because many x/z levels\n")}
  }
  
  dodgeVal <- 0.6
  colAlpha <- 1
  jitterSize <- 0.02 # used to be 0.05
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  
  cat("Overall condition means (without first aggregating per subject):\n")
  print(tapply(data[, yVar], interaction(data[, zVar], data[, xVar]), mean, na.rm = T))
  
  cat("Create new variables x, y, z, subject based on inputs\n")
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$z <- data[, zVar]
  data$subject <- data[, subVar]
  
  ## Exclude NAs:
  nRow1 <- nrow(data)
  data <- droplevels(subset(data, !(is.na(x)) & !(is.na(y)) & !(is.na(z))))
  nRow2 <- nrow(data)
  cat(paste0("Excluded ", nRow1 - nRow2, " rows due to NAs\n"))
  
  ## Colours:
  if(is.null(selCol)){
    selCol <- retrieve_colour(zVar)
    selCol <- rep(selCol, length.out = length(unique(data[, zVar])))
  }
  if(length(selCol) != length(unique(data[, zVar]))){
    stop(paste0("Length selCol = ", length(selCol), " while number levels zVar = ", length(unique(data[, zVar])), ", do not match"))
  }
  condCol <- rep(selCol, length(unique(data[, xVar])))
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  cat("Aggregate data per subject\n")
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/condition1/condition2 in one line, variables subject, x, y, z
  
  ## Add condition variable:
  cat("Create condition variable\n")
  nZlevel <- length(unique(data$z))
  posScale <- 0.05 * (nZlevel + 1)
  aggrData$cond <- as.numeric(aggrData$x)*nZlevel - nZlevel + as.numeric(aggrData$z) # necessary for proper axis positions
  # aggrData$cond <- 1 + as.numeric(aggrData$x)*2 + as.numeric(aggrData$z)
  nCond <- length(unique(aggrData$cond))
  if (length(condCol) < nCond){condCol <- rep(condCol, length.out = nCond)}
  
  ## Add jittered x-axis for points:
  cat("Add jitter for points\n")
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = jitterSize) # pure jitter .05
  if (nZlevel == 2){
    aggrData$xpos <- as.numeric(aggrData$x) - min(as.numeric(aggrData$x)) + 1 + (as.numeric(aggrData$z) - 1.5) * 2 * posScale # convert to [1 2], to [-0.5 0.5], * 2 so [-1 1], scale by 0.15
  } else if (nZlevel == 3) {
    zMid <- round(mean(as.numeric(data$z)))
    aggrData$xpos <- as.numeric(aggrData$x) - min(as.numeric(aggrData$x)) + 1 + ((as.numeric(aggrData$z) - zMid)) * posScale  # demean, scale by 0.20
  } else {
    zMid <- round(mean(as.numeric(data$z)))
    zScale <- ceiling(max(as.numeric(data$z))) - zMid
    aggrData$xpos <- as.numeric(aggrData$x) - min(as.numeric(aggrData$x)) + 1 + ((as.numeric(aggrData$z) - zMid)) * zScale * posScale  # demean, bring min/max to 1, scale by 0.20
    warning(paste0("Not yet implement for z variable with ", nZlevel, " levels\n"))
  }
  aggrData$xj <- aggrData$xpos + aggrData$j # add jitter to xpos
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    cat("Automatically determine y-axis limits based on per-subject-per-condition means\n")
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  cat("Aggregate data across subjects\n")
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x", "z"))
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  ## Additions:
  # prefix <- NULL; suffix <- NULL
  if(is.null(prefix)){prefix <- ""} else {prefix <- paste0(prefix, "_")}
  if(is.null(suffix)){suffix <- ""} else {suffix <- paste0("_", suffix)}
  
  ## Name:
  plotName <- paste0("custombarplot2_", prefix, yVar, "_", xVar, "_", zVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  if (isBeeswarm){plotName <- paste0(plotName, "_beeswarm")} 
  plotName <- paste0(plotName, suffix)
  cat(paste0("Start plot ", plotName, "\n"))
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  ## Start plot:
  p <- ggplot(summary_d, aes(x, y, fill = z))
  
  ## Bars of means:
  cat("Add bars\n")
  p <- p + stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
                        lwd = LWD, color = "black")
  
  ## Error bars:
  cat("Add error bars\n")
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = 0.2,
                         lwd = LWD, color = "black", alpha = 1)
  
  # Individual data points:
  if (isPoint){
    cat("Start adding per-subject points \n")
    for(iCond in 1:nCond){ # add separately per condition
      p <- p + geom_point(data = aggrData[aggrData$cond == iCond, ],
                          aes(x = xj), # position = "dodge",
                          shape = 21, size = 2, stroke = 1.2, color = "black", fill = condCol[iCond],
                          alpha = 0.5) # colAlpha)
    }
  }
  
  ## Connect colored dots:
  if (isConnect){
    cat("Start connecting per-subject points \n")
    
    subVec <- sort(unique(aggrData$subject))
    nSub <- length(subVec)
    
    for(iSub in 1:nSub){ # iSub <- 1
      subData <- subset(aggrData, subject == subVec[iSub])
      p <- p + geom_path(data = subData, aes(x = xj, y = y, group = 1), color = 'grey40', # color = 'grey70'
                         alpha = 0.50, size = dotSize) # consider alpha = 0.80
    }
  }  
  
  ## Beeswarm style plots:
  if (isBeeswarm){
    cat("Start adding beeswarm \n")
    for(iCond in 1:nCond){ # add separately per condition
      p <- p + geom_beeswarm(data = aggrData[aggrData$cond == iCond, ],
                             aes(x = xpos), # position = "dodge",
                             # priority = "ascending",
                             shape = 21, size = 2, stroke = 1.2, color = "black", fill = condCol[iCond],
                             alpha = 0.5) # colAlpha)
    }
  }
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  # Settings:
  if (yLim[1] == 0 & yLim[2] == 1){
    cat("Add y-axis ticks for 0, 0.5, 1\n")
    # p <- p + scale_y_break(c(0, 0.5, 1))
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim = yLim)}
  # if(!(is.null(yLim))){p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1], yLim[-1], (yLim[-1] - yLim[1])/2))}
  
  # Add theme, font sizes:
  cat("Add axis labels, colors, theme, font sizes\n")
  require(ggthemes)
  p <- p + labs(x = xLab, y = yLab, fill = zLab) +
    scale_fill_manual(values = selCol, limits = levels(summary_d$z)) + 
    theme_classic() + 
    theme(axis.text = element_text(size = FTS),
          axis.title = element_text(size = FTS), 
          plot.title = element_text(size = FTS, hjust = 0.5), 
          axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
  
  ## Add legend:
  if (addLegend){
    p <- p + theme(
      legend.text = element_text(size = FTS),
      legend.title = element_text(size = FTS)
    )
  } else {
    p <- p + theme(
      legend.title = element_blank(), legend.position = "none"
    )
  }
  
  print(p)
  if(savePNG | saveEPS){dev.off(); print(p)}
  cat("Finished :-)\n")
  return(p)
}

# ================================================================================================================================================ #
#### Barplot 3 IVs: Aggregate per condition per subject, plot (3 IVs, 1 on x-axis, 1 as color, 1 as facets): ####

custom_barplot3 <- function(data, yVar, xVar, zVar, splitVar = NULL, subVar = "subject_n", 
                            yLab = NULL, xLab = NULL, zLab = NULL, main = NULL,
                            xLevels = NULL, zLevels = NULL, splitLevels = NULL,
                            selCol = NULL, isPoint = T, yLim = NULL, savePNG = T, saveEPS = F){
  #' Make bar plots with 3 IVs: x-axis and color and facetwrap.
  #' Can add points with geom_point, 
  #' but not beeswarm plots because no position argument (hence no dodge) and manual x-position not compatible with facet_wrap.
  #' Note: In order to get the order of bars (from left to right) correctly, all factors are recoded into values from
  #' (nLevels - 1) to 0 in descending order. Later, factor levels are added again in the correct order.
  #' Check in print out if relabeling is done correctly!
  #' @param data data frame, trial-by-trial data.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param zVar string, name of variable that determines bar coloring. Needs to be a factor.
  #' @param splitVar string, name of variable by which to split plot (facetwrap, optional).
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: retrieve appropriate name with substitute_label()).
  #' @param yLab string, label for y-axis (default: retrieve appropriate name with substitute_label()).
  #' @param zLab string, label for color legend (default: retrieve appropriate name with substitute_label()).
  #' @param main string, title of plot (optional).
  #' @param xLevels string, level names for x-axis (default: retrieve from xVar in alphabetical order).
  #' @param zLevels string, level names for x-axis (default: retrieve from zVar in alphabetical order).
  #' @param splitLevels string, level names for x-axis (default: retrieve from splitVar in alphabetical order).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @return creates (and saves) plot.  #' Make bar plot per subject
  #' @param data data frame, with variables \code{variables}
  
  # yLab = NULL; xLab = NULL; zLab = NULL; main = NULL;
  # xLevels = NULL; zLevels = NULL; splitLevels = NULL;
  # selCol = NULL; isPoint = F; yLim = NULL; savePNG = T; saveEPS = F
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  require(ggbeeswarm) # for ggbeeswarm
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## General plotting settings:
  
  LWD <- retrieve_plot_defaults("LWD") # 1.3 # 1.5
  if (savePNG | saveEPS){FTS <- retrieve_plot_defaults("FTS")} else {FTS <- 15} # 30
  dodgeVal <- 0.6
  barWidth <- 0.15
  colAlpha <- 0.6
  
  SEweight <- 1.96
  
  isPrint <- T # print data sets to check proper recoding
  
  ## Set default y-axis limits:
  if (is.null(yLim)){
    yLim <- c(0, 1)
  }
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  ## Input variables:
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(zVar %in% names(data))){stop("zVar not found in data")}
  if(!is.null(splitVar) & !(subVar %in% names(data))){stop("splitVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  
  ## Retrieve axis labels:
  if(is.null(xLab)){xLab <- substitute_label(xVar)}
  if(is.null(yLab)){yLab <- substitute_label(yVar)}
  if(is.null(zLab)){zLab <- substitute_label(zVar)}
  zLab <- gsub(" ", " \n", zLab) # add newline to any z-label
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  
  cat("Create new numerical variables x, y, z, subject based on inputs\n")
  data$y <- data[, yVar]
  data$x <- as.numeric(data[, xVar]) # convert to numeric
  data$z <- as.numeric(data[, zVar]) # convert to numeric
  if (!is.null(splitVar)){data$split <- as.numeric(data[, splitVar])} # convert to numeric
  data$subject <- data[, subVar]
  
  ## Exclude NAs:
  nRow1 <- nrow(data)
  data <- droplevels(subset(data, !(is.na(x)) & !(is.na(y)) & !(is.na(z)) & !(is.na(split))))
  nRow2 <- nrow(data)
  cat(paste0("Excluded ", nRow1 - nRow2, " rows due to NAs\n"))
  
  ## Determine level names if not set as input:
  if(is.null(xLevels)){xLevels <- sort(as.character(unique(data[, xVar])))}
  if(is.null(zLevels)){zLevels <- sort(as.character(unique(data[, zVar])))}
  if (!all(levels(data[, zVar]) == sort(levels(data[, zVar])))){zLevels <- rev(zLevels); cat("Factor levels of zVar not in alphabetical order, invert\n")} # invert levels of z if factor levels not in alphabetical order
  if(is.null(splitLevels)){
    splitLevels <- sort(unique(data[, splitVar]))
    # if (!all(levels(data[, splitVar]) == sort(levels(data[, splitVar])))){splitLevels <- rev(splitLevels); cat("Factor levels of splitVar not in alphabetical order, invert\n")} # invert levels of split if factor levels not in alphabetical order
  }
  
  ## Colours (must be determined after NA exclusion):
  if(is.null(selCol)){
    selCol <- retrieve_colour(zVar)
    selCol <- rep(selCol, length.out = length(unique(data[, zVar]))) # copy as often as necessary
  }
  if(length(selCol) != length(unique(data[, zVar]))){
    stop(paste0("Length selCol = ", length(selCol), " while number levels zVar = ", length(unique(data[, zVar])), ", do not match"))
  }
  condCol <- rep(selCol, length(unique(data[,xVar])))
  
  # -------------------------------------------------------------------------- #
  ### Recode to 0 until (nLevels - 1) for combining into condition variable:
  
  ## X-axis levels:
  xLevelVec <- sort(unique(data$x)) # levels
  nXlevels <- length(xLevelVec) # number of levels
  if (isPoint & nXlevels > 3){stop("Points not implemented for x variable with > 3 levels")}
  
  ## Z-axis levels:
  zLevelVec <- sort(unique(data$z)) # levels
  nZlevels <- length(zLevelVec) # number of levels
  if (isPoint & nZlevels > 2){stop("Points not implemented for z variable with > 2 levels")}
  
  ## Recode variables to go from (nLevels - 1) till 0 in descending order:
  data$x <- nXlevels - data$x # recode to (nXlevels - 1) to 0 in descending order
  data$z <- nZlevels - data$z # recode to (nXlevels - 1) to 0 in descending order
  xLevelVec <- rev(sort(unique(data$x))) # update, descending order
  zLevelVec <- rev(sort(unique(data$z))) # update, descending order
  
  if (!is.null(splitVar)){
    sLevelVec <- sort(unique(data$split))
    nSlevels <- length(sLevelVec)
    if (isPoint & nSlevels > 3){stop("Points not implemented for split variable with > 2 levels")}
    data$split <- max(data$split) - data$split # recode to (nXlevels - 1) to 0 in reverse order
    sLevelVec <- rev(sort(unique(data$split))) # update, descending order
  } # to 0 - 1
  
  ## Print recoding of factor levels into descending numbers to console: 
  if(isPrint){print(table(data[, xVar], data$x))} # 1-N becomes (N-1)-0
  if(isPrint){print(table(data[, zVar], data$z))} # 1-N becomes (N-1)-0
  if(isPrint & !is.null(splitVar)){print(table(data[, splitVar], data$split))} # # 1-N becomes (N-1)-0
  
  # -------------------------------------------------------------------------- #
  ## Combine into single condition variable:
  
  if (!is.null(splitVar)){ # if splitVar: 8 conditions
    data$condition <- 1 + data$split*nXlevels*nZlevels + data$x*nZlevels + data$z
  } else { # No splitVar: 4 conditions
    data$condition <- 1 + data$x*nZlevels + data$z
  }
  nCond <- length(unique(data$condition))
  stopifnot(min(data$condition)  == 1)
  stopifnot(max(data$condition)  == nCond)
  
  # sort(unique(data$condition)) # from 1 to nCond
  # table(data$condition, data$split) # slowest factor, everything nXlevels*nZlevels times, from 0-nSlevels in ascending order
  # table(data$condition, data[, splitVar]) # slowest factor, everything nXlevels*nZlevels times, from nSlevels-0 in descending order
  # table(data$condition, data$x) # middle factor, everything nZlevels times, from 0-nXlevels in ascending order
  # table(data$condition, data[, xVar]) # middle factor, everything nZlevels times, from nXlevels-0 in descending order
  # table(data$condition, data$z) # fastest factor, odd and even, from 0-nZlevels in ascending order
  # table(data$condition, data[, zVar]) # fastest factor, odd and even, from nZlevels-0 in descending order
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  
  cat("Aggregate data per subject\n")
  aggrData <- ddply(data, .(subject, condition), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  
  ## Recover original variables:
  aggrData$x <- (ceiling(aggrData$condition/nZlevels) - 1) %% nXlevels # middle factor, everything nZlevels times 
  aggrData$z <- (aggrData$condition - 1) %% nZlevels # fastest factor, odd or even
  if (!is.null(splitVar)){
    aggrData$split <- ceiling(aggrData$condition/(nXlevels*nZlevels) - 1)
  }
  
  # print(head(aggrData, n = nCond)) # needs to match condition meaning in data object above
  # table(aggrData$condition, aggrData$split) # slowest factor, everything nXlevels*nZlevels times, from 0-nSlevels in ascending order
  # table(aggrData$condition, aggrData$x) # middle factor, everything nZlevels times, from 0-nXlevels in ascending order
  # table(aggrData$condition, aggrData$z) # fastest factor, odd and even, from 0-nZlevels in ascending order
  
  ### Create factors:
  ## Reverse above inversion: xLevelVec is inverted, but xLevels the right way around
  aggrData$x_f <- factor(aggrData$x, levels = xLevelVec, labels = xLevels) # assign factor levels to numerics (in descending order) 
  aggrData$z_f <- factor(aggrData$z, levels = zLevelVec, labels = zLevels) # assign factor levels to numerics (in descending order)
  if (!is.null(splitVar)){aggrData$split_f <- factor(aggrData$split, levels = sLevelVec, labels = splitLevels)} # assign factor levels to numerics (in descending order)
  
  cat(paste0("Assume ", paste0(xLevelVec, collapse = ", "), " corresponds to ", paste0(xLevels, collapse = ", "), "\n"))
  cat(paste0("Assume ", paste0(zLevelVec, collapse = ", "), " corresponds to ", paste0(zLevels, collapse = ", "), "\n"))
  if (!is.null(splitVar)){cat(paste0("Assume ", paste0(sLevelVec, collapse = ", "), " corresponds to ", paste0(splitLevels, collapse = ", "), "\n"))}
  
  if(isPrint){print(head(aggrData, n = nCond))} # needs to match condition meaning in data object above
  
  # -------------------------------------------------------------------------- #
  ## Determine y limits if not given:
  
  if(is.null(yLim)){
    cat("Automatically y-axis limits based on per-subject-per-condition means\n")
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  cat("Aggregate data across subjects\n")
  
  d <- summarySEwithin(aggrData, measurevar = "y", withinvar = "condition", idvar = "subject", na.rm = T)
  d$condition <- as.numeric(as.factor(d$condition)) # condition back to numeric
  
  ## Recover independent variables from condition:
  d$x <- ceiling(d$condition/nZlevels - 1) %% nXlevels
  d$z <- (d$condition - 1) %% nZlevels
  if (!is.null(splitVar)){
    d$split <- ceiling(d$condition/(nXlevels*nZlevels) - 1)
  }
  
  # print(d)
  
  ## Create factors:
  d$x_f <- factor(d$x, levels = xLevelVec, labels = xLevels)
  d$z_f <- factor(d$z, levels = zLevelVec, labels = zLevels)
  if (!is.null(splitVar)){d$split_f <- factor(d$split, levels = sLevelVec, labels = splitLevels)}
  
  if(isPrint){print(d)}
  
  ## Check if y +/- se within ylim:
  if (any(d$y - d$se < yLim[1])){warning("Lower error bars will exceed y-axis limit")}
  if (any(d$y + d$se > yLim[2])){warning("Upper error bars will exceed y-axis limit")}
  
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ## Start plot:
  
  ## Name:
  plotName <- paste0("custombarplot3_",  yVar, "~", xVar, "_", zVar, "_")
  if (!(is.null(splitVar))){plotName <- paste0(plotName, splitVar)}
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  cat(paste0("Start plot ", plotName, "\n"))
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plot, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plot, plotName, ".png"), width = 480, height = 480)}
  
  ## Initialize ggplot object:
  p <- ggplot(d, aes(x = x_f, y = y, fill = z_f))
  
  ## Bars of means:
  cat("Add bars\n")
  p <- p + geom_bar(position = "dodge", stat = "summary", fun = "identity",
                    color = "black", width = dodgeVal, lwd = LWD)
  
  ## Error bars:
  cat("Add error bars\n")
  p <- p + geom_errorbar(data = d,
                         aes(x = x_f, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = barWidth,
                         lwd = LWD, color = "black", alpha = 1)
  
  ## Individual data points:
  if (isPoint){
    cat("Start adding per-subject points \n")
    p <- p + geom_point(data = aggrData,
                        position = position_dodge(width = dodgeVal),
                        shape = 21, size = 2, stroke = 1.2, color = "black",
                        alpha = 0.5) 
  }
  
  ## Facet wrap:
  if (!is.null(splitVar)){
    cat("Start adding facet_wrap\n")
    p <- p + facet_wrap(vars(split_f))
  }
  
  ## Y axis limits:
  cat("Start adding y-axis ticks\n")
  p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1], yLim[-1], (yLim[-1] - yLim[1])/2)) 
  
  ## Add labels:
  cat("Start labels\n")
  p <- p + labs(x = xLab, fill = zLab,
                y = yLab)
  
  ## Add color:
  cat("Start colors for fill\n")
  p <- p + scale_fill_manual(values = rep(selCol, 4), limits = levels(d$z_f))
  
  ## Add theme:
  cat("Start theme\n")
  p <- p + theme_classic()
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  ## Font sizes:
  cat("Start line width and font size \n")
  p <- p + theme(axis.line = element_line(colour = 'black'), # linewidth = LWD),
                 axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5), # center title 
                 legend.title = element_text(size = FTS),
                 strip.text.x = element_text(size = FTS), # facetwrap FTS
                 legend.text = element_text(size = FTS))
  print(p)
  if (savePNG | saveEPS){ dev.off(); print(p); Sys.sleep(1)}
  cat("Finished :-)\n")
  
  return(p)
  
} # end of function

# ================================================================================================================================================ #
#### Lineplot 1 IV: Aggregate per time point per condition per subject, plot (1 IV for any condition, time on x-axis): ####

custom_lineplot <- function(data, xVar = "counter", yVar = "response_cleaned", zVar = "condition_f", subVar = "subject_n", 
                            xLab = "Time (trial number)", yLab = "p(Go)", main = "",
                            selCol = c("#009933", "#CC0000", "#009933", "#CC0000"), selLineType = c(1, 1, 2, 2),
                            SEweight = 1, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make line plot with group-level and individual lines.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. Variable needs to be numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Variable needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (optional).
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: c("#009933", "#CC0000", "#009933", "#CC0000")).
  #' @param selLineType vector of numerics, line types to use (default: c(1, 1, 2, 2))
  #' @param SEweight scalar, weight to use for error shades (how many times SE; default: 1).
  #' @param yLim vector of two numbers, y-axis (default: NULL).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  LWD <- retrieve_plot_defaults("LWD") * 2 # 3 # axes of plot
  CEX <- 1.5 # axes ticks and labels
  lineWidth <- retrieve_plot_defaults("LWD") * 2 # 3
  FTS <- retrieve_plot_defaults("FTS")
  dodgeVal <- retrieve_plot_defaults("dodgeVal")
  colAlpha <- 1
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  data$x <- data[,xVar]
  data$y <- data[,yVar]
  data$z <- data[,zVar]
  data$subject <- data[,subVar]
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data:
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/x condition/z condition in one line, variables subject, x, y, z
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x", "z"))
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # -------------------------------------------------------------------------- #
  ## Data dimensions:
  xVec <- unique(sort(as.numeric(summary_d$x)))
  xMax <- max(xVec)
  condNames <- unique(summary_d$z)
  nCond <- length(unique(summary_d$z))
  
  # -------------------------------------------------------------------------- #
  ## Plot name:
  
  plotName <- paste0("lineplot_",yVar, "_",xVar, "_",zVar)
  
  # -------------------------------------------------------------------------- #
  # Saving:
  
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plot, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plot, plotName, ".png"), width = 480, height = 480)}
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  par(mar = c(5.1, 5.1, 4.1, 2.1)) # bottom, left, top, right
  
  # dev.off()
  for (iCond in 1:nCond){ # iCond <- 1
    condName <- condNames[iCond] # name of condition
    yVec <- summary_d$y[summary_d$z == condName] # y-variable
    seVec <- summary_d$se[summary_d$z == condName] # se variable
    plot(xVec, yVec, type = "l", 
         col = selCol[iCond], lty = selLineType[iCond], axes = F,
         lwd = LWD, cex.lab=CEX, cex.axis=CEX, cex.main=CEX,
         xlab = xLab, ylab = yLab, main = main,
         xlim = c(0, xMax), ylim = yLim)
    axis(side = 1, lwd = LWD, cex.axis = CEX, at = seq(0,xMax, 5), line = 0)
    axis(side = 2, lwd = LWD, cex.axis = CEX, at = c(0, 0.5, 1))
    polygon(c(xVec,rev(xVec)),
            c(yVec-SEweight*seVec,rev(yVec+SEweight*seVec)),col = alpha(selCol[iCond],0.2), border = F)
    par(new = TRUE)
    
  }
  
  # Add legend:
  legend("top", legend=condNames,
         col=selCol, lty = selLineType, border = 0, lwd = LWD, cex = CEX, horiz=TRUE, bty = "n")
  
  if(savePNG | saveEPS){dev.off()}
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right
  
}

# ================================================================================================================================================ #
#### Lineplot 1 IV with ggplot: Aggregate per time point per condition per subject, plot (1 IV for any condition, time on x-axis): ####

custom_lineplot_gg <- function(data, xVar, yVar, zVar, subVar = "subject_n", 
                               xLab = NULL, yLab = NULL, main = NULL,
                               selCol = NULL, selLineType = c(1), breakVec = NULL,
                               SEweight = 1, yLim = NULL, addLegend = F, savePNG = F, saveEPS = F){
  #' Make line plot with group-level lines plus shades in ggplot.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. Variable needs to be numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Variable needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: retrieve appropriate name with substitute_label()).
  #' @param yLab string, label for y-axis (default: retrieve appropriate name with substitute_label()).
  #' @param main string, overall plot label (optional).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param selLineType vector of numerics, line types to use (default: c(1, 1, 2, 2))
  #' @param SEweight scalar, weight to use for error shades (how many times SE; default: 1).
  #' @param yLim vector of two numbers, y-axis (default: NULL).
  #' @param addLegend Boolean, add legend at the top or not (default: FALSE).
  #' @param savePNG Boolean, save as .png file (default: FALSE).
  #' @param saveEPS Boolean, save as .eps file (default: FALSE).
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  ## Input variables:
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  
  ## Axis labels:
  if(is.null(xLab)){xLab <- substitute_label(xVar)}
  if(is.null(yLab)){yLab <- substitute_label(yVar)}
  # if(is.null(zLab)){zLab <- substitute_label(zVar)}
  
  ## Colours:
  if(is.null(selCol)){selCol <- retrieve_colour(zVar)}
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  LWD <- retrieve_plot_defaults("LWD") # 1.3
  FTS <- retrieve_plot_defaults("FTS")
  colAlpha <- 1
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$z <- data[, zVar]
  data$subject <- data[, subVar]
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per conditions x/y:
  
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/x condition/z condition in one line, variables subject, x, y, z
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data across subjects with Rmisc:
  
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x", "z"))
  summary_d$x <- as.numeric(summary_d$x) # back to numeric to get continuous x-axis
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # Data dimensions:
  xVec <- unique(sort(as.numeric(summary_d$x)))
  condNames <- unique(summary_d$z)
  nCond <- length(unique(summary_d$z))
  
  if (length(selLineType) == 1){selLineType <- rep(selLineType, nCond)}
  
  # -------------------------------------------------------------------------- #
  ## Start plot with ggplot:
  
  ## Plot name:
  plotName <- paste0("lineplot_gg_", yVar, "_", xVar, "_", zVar)
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plot, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plot, plotName, ".png"), width = 480, height = 480)}
  
  ## Start plot:
  # par(mar = c(5.1, 5.1, 4.1, 2.1)) # bottom, left, top, right
  p <- ggplot(data = summary_d, aes(x = x, y = y, col = z, fill = z))
  
  # p + geom_path() # just to get legend going
  
  ## Loop over conditions, make shade and line on top:
  for (iCond in 1:nCond){ # iCond <- 1
    
    ## Select data, set upper and lower shade limits:
    condData <- subset(summary_d, z == condNames[iCond]) # select data for this condition
    condData$ymin <- condData$y - SEweight * condData$se # lower edge of shade
    condData$ymax <- condData$y + SEweight * condData$se # upper edge of shade
    
    ## Shade:
    p <- p + geom_ribbon(data = condData, aes(x = x, y = y, ymin = ymin, ymax = ymax, group = 1), 
                         # col = NA, # remove outer border of shades
                         linetype = 0, # remove outer border of shades
                         # fill = selCol[iCond], # comment out for legend in colours (not grey)
                         # show.legend = T, 
                         alpha = 0.2, lwd = 0)
    
    ## Line:
    p <- p + geom_path(data = condData, aes(x = x, y = y, group = 1),
                       # col = selCol[iCond], 
                       linetype = selLineType[iCond], size = LWD) 
  }
  
  # Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  ## X-axis:
  xMin <- min(xVec)
  xMax <- max(xVec)
  if (is.null(breakVec)){
    xRange <- xMax - xMin
    xStep <- ceiling(xRange/5)
    breakVec <- seq(xMin, xMax, xStep)
  }
  p <- p + scale_x_continuous(limits = c(xMin, xMax), breaks = breakVec)
  
  ## Y-axis:
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim=yLim)}
  
  ## Add labels:
  p <- p + labs(x = xLab, y = yLab, fill = condNames, col = condNames)
  
  ## Add line colors:
  p <- p + scale_fill_manual(values = selCol, limits = levels(data$z)) # set limits for correct order
  p <- p + scale_color_manual(values = selCol, limits = levels(data$z), guide = "none") # set limits for correct order
  
  # Add theme, font sizes:
  require(ggthemes)
  p <- p + theme_classic() + 
    theme(axis.text = element_text(size = FTS),
          axis.title = element_text(size = FTS), 
          plot.title = element_text(size = FTS, hjust = 0.5))
  
  ## Add legend:
  if (addLegend){
    p <- p + theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = element_blank(),
      # legend.title = element_text(size = FTS),
      legend.text = element_text(size = FTS/2)
    )
  } else {
    p <- p + theme(
      legend.title = element_blank(), legend.position = "none"
    )
  }
  
  print(p)
  if(savePNG | saveEPS){dev.off()}
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right
  
  return(p)
  
}

# ===============================================================================================#
#### Plot density split by zVar: ####

customplot_density2 <- function(data, xVar, zVar, 
                                xLab = NULL, zLab = NULL, xLim = NULL, main = NULL,
                                selCol = NULL, addLegend = F, isPNG = T){
  #' Plot density of xVar split by yVar.
  #' @param data data frame, trial-by-trial data.
  #' Param xVar scalar string, variable for which to compute density.
  #' @param zVar scalar string, name of variable to split by (differently colored lines). Variable needs to be a factor.
  #' @param xLab scalar string, x-axis label (optional).
  #' @param zLab scalar string, color label (optional).
  #' @param main string, overall plot label (optional).
  #' @param xLim vector of two numbers, x-axis limits (optional).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param addLegend Boolean, add legend at the top or not (default: FALSE).
  #' @param isPNG Boolean, save as .png file (default: TRUE).
  #' @return creates (and saves) plot.
  
  require(ggplot2)
  cat(paste0("Plot density of ", xVar, " split by ", zVar, "\n"))
  
  if(!is.numeric(data[, xVar])){stop("xVar must be numeric")}
  if(is.numeric(data[, zVar])){stop("zVar must be numeric")}
  
  # -------------------------------------------------------------------------- #
  ### Retrieve settings:
  
  FTS <- retrieve_plot_defaults("FTS")
  LWD <- retrieve_plot_defaults("LWD")
  if (is.null(xLab)){xLab <- substitute_label(xVar)}
  if (is.null(zLab)){zLab <- substitute_label(zVar)}
  if (is.null(selCol)){
    selCol <- retrieve_colour(zVar)
  } else {
    if(length(selCol) != length(levels(data[, zVar]))){"selCol and level(data$zVar) of different lengths"}
  }
  
  if (is.null(xLim)){
    if (grepl("RT_n", xVar)){
      xLim <- c(0, 1.3)
      # } else if (grepl("dilation_n", xVar)){
      #     xLim <- c(0, 80)
    } else {
      xLim <- c(min(data[, xVar], na.rm = T), max(data[, xVar], na.rm = T))
    }
  }
  
  ## Copy over data:
  data$x <- data[, xVar]
  data$z <- data[, zVar]
  
  ## Exclude NAs:
  data <- data[which(!is.na(data$x) & !is.na(data$z)), ]  
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  p <- ggplot(data = data, aes(x = x, fill = z)) # initialize
  p <- p + geom_density(color = "black", alpha = 0.7, trim = FALSE, lwd = LWD)
  
  ## Fill colour:
  p <- p + scale_fill_manual(values = selCol)
  
  ## Limits:
  p <- p + coord_cartesian(xlim = xLim) 
  
  # -------------------------------------------------------------------------- #
  ### Labels:
  
  ## Labels:
  p <- p + labs(x = xLab, y = "Density", fill = zLab)
  
  ## Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  # -------------------------------------------------------------------------- #
  ### Add theme, font sizes:
  
  ## Theme:
  p <- p + theme_classic() # theme
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = FTS),
                 axis.title = element_text(size = FTS), 
                 plot.title = element_text(size = FTS, hjust = 0.5)) # center title 
  
  ## Add legend:
  if (addLegend){
    p <- p + theme(
      # legend.title = element_blank(),
      legend.title = element_text(size = FTS),
      legend.text = element_text(size = FTS)
    )
  } else {
    p <- p + theme(
      legend.title = element_blank(), legend.position = "none"
    )
  }
  
  # -------------------------------------------------------------------------- #
  ### Save:
  if (isPNG){
    ## Save:  
    figName <- paste0("density_", xVar, "~", zVar)
    if(addLegend){figName <- paste0(figName, "_addLegend")}
    png(paste0(dirs$plotDir, figName,  ".png"), width = 480, height = 480)
    print(p)
    dev.off()
  }
  
  # -------------------------------------------------------------------------- #
  ### Return:
  
  print(p)
  return(p)
  cat("Finished! :-)\n")
  
}

# ================================================================================================================================================ #
#### Plot correlation (scatterplot and regression line) between two variables: ####

plot_correlation <- function(data, xVar, yVar, 
                             isSubLabel = F, subVar = "subject_n", 
                             xLab = NULL, yLab = NULL, main = NULL, printCor = T,
                             xLim = NULL, yLim = NULL, FTS = NULL, isSave = F, useLaTeX = F){
  #' Plot correlation (scatterplot and regression line) between two variables in data frame with ggplot.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. Variable needs to be numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param isSubLabel Boolean, whether to add subject IDs as labels next to dots (T) or not (F) (default: F).
  #' @param subVar string, name of subject identifier used for labeling points. Variable needs to be numeric.
  #' @param xLab string, label for x-axis (default: "X").
  #' @param yLab string, label for y-axis (default: "Y").
  #' @param main string, title (default: none).
  #' @param printCor Boolean, if no main provided, print correlation or not (default: TRUE).
  #' @param xLim vector of two numerics, x-axis limits (optional).
  #' @param yLim vector of two numerics, y-axis limits (optional).
  #' @param FTS numeric, font size of axes (optional).
  #' @param isSave Boolean, save as .png file (T) or not (F; default: F).
  #' @return creates (and saves) plot.
  
  # isSubLabel = F; subVar = "subject_n"; xLab = NULL; yLab = NULL; main = NULL; xLim = NULL; yLim = NULL; FTS = NULL; isSave = F
  
  require(ggplot2)
  require(latex2exp)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ### Fixed settings:
  if (is.null(FTS)){
    FTS <- 20 
  }
  LWD <- 1.3
  MKS <- 3
  subFTS <- 6 # 8
  nRound <- 3 # for print xMin/ xMax/ yMin/ yMax to console
  # margin <- c(1.0, 1.5, 1.0, 1.0, "cm") # top, right, bottom, left
  
  nData <- sum(complete.cases(data[, c(xVar, yVar)]))
  cat(paste0("Plot correlation between ", xVar, " and ", yVar, " given ", nData, " data points\n"))
  
  ## Axis labels:
  if (is.null(xLab)){xLab <- xVar}
  if (is.null(yLab)){yLab <- yVar}
  
  ## Print correlation as header:
  if (is.null(main) & printCor){
    r <- as.numeric(cor(data[xVar], data[yVar], use = "complete.obs"))
    if (r < 0){nPad <- 6} else {nPad <- 5}
    rText <- str_pad(round(r, nRound), nPad, "right", "0")
    df <- nData - 2
    t <- r/sqrt((1 - r^2) / df)
    pVal <- pt(abs(t), df, lower.tail = F) * 2
    if (pVal < 0.001){
      pValText <- "p < .001"
    } else {
      pValText <- paste0("p = .", str_pad(round(pVal, nRound)*1000, nRound, "left", "0"))
    }
    main <- paste0("r(", df, ") = ", rText, ", ", pValText)
    cat(paste0(main, "\n"))
  }
  
  # -------------------------------------------------------------------------- #
  ### Name for saving:
  
  if (isSave){
    plotName <- paste0("ggplot_correlation_", xVar, "_", yVar, "_", FTS, ".png")
    cat(paste0("Save plot under ", plotName, "\n"))
    png(paste0(dirs$plot, plotName),
        width = 480, height = 480) # FTS <- 17
    # tiff(paste0(dirs$plot, "ggplot_correlation_", xVar, "_", yVar, "_", FTS, ".tiff"),
    #     width = 480, height = 320) # FTS <- 17
  }
  
  # -------------------------------------------------------------------------- #
  ### Start ggplot:
  
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  
  if (isSubLabel){
    data$subject <- data[, subVar]
    p <- ggplot(data, aes(x = x, y = y, label = subject)) + 
      geom_text(hjust = 0, vjust = 0, size = subFTS) # hjust = -0.5
  } else {
    p <- ggplot(data, aes(x = x, y = y))
  }
  
  ### Points for scatter plot:
  p <- p + geom_point(shape = 19, size = MKS, stroke = 0.8, fill = "white")
  
  ### Regression line:
  p <- p + geom_smooth(method = lm, se = T, col = "red", linewidth = 2) # regression line 
  
  ### Horizontal and vertical lines:
  p <- p + geom_vline(xintercept = 0, linetype = 2, color = "black", size = 1) # Middle line at 0
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "black", size = 1) # Middle line at 0
  
  ### Labels:
  if (useLaTeX){
    # p <- p + labs(x = TeX(paste0("$", xLab, "$")), y = TeX(paste0("$", yLab, "$")), title = "")
    p <- p + labs(x = TeX(paste0(xLab)), y = TeX(paste0(yLab)), title = "")
  } else {
    p <- p + xlab(xLab) + ylab(yLab) # Labels
  }

  ### Title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)
  }
  
  ### Theme:
  p <- p + theme_classic() + # theme
    theme(plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"), # top, right, bottom, left
          axis.text = element_text(size = FTS, color = "black"),
          axis.title = element_text(size = FTS, color = "black"), 
          plot.title = element_text(size = FTS, color = "black", hjust = 0.5), # center title 
          legend.text = element_text(size = FTS, color = "black"),
          axis.line = element_line(colour = 'black', linewidth = LWD))
  
  # -------------------------------------------------------------------------- #
  ### Axis limits:
  
  ## Retrieve:
  xMin <- min(data[, xVar])
  xMax <- max(data[, xVar])
  yMin <- min(data[, yVar])
  yMax <- max(data[, yVar])
  cat(paste0("X ranges from ", round(xMin, nRound), " to ", round(xMax, nRound), "\n"))
  cat(paste0("Y ranges from ", round(yMin, nRound), " to ", round(yMax, nRound), "\n"))
  
  ## Erode:
  # if(xMin > 0){xMin <- xMin * 0.90} else (xMin <- xMin * 1.10)
  # if(xMax > 0){xMax <- xMax * 1.10} else (xMax <- xMax * 0.90)
  # if(yMin > 0){yMin <- yMin * 0.90} else (yMin <- yMin * 1.10)
  # if(yMax > 0){yMax <- yMax * 1.10} else (yMax <- yMax * 0.90)
  
  ## Set limits:
  if(is.null(xLim)){xLim <- c(xMin, xMax); cat("Set xLim automatically\n")} # define x-axis limits
  if(is.null(yLim)){yLim <- c(yMin, yMax); cat("Set yLim automatically\n")} # define x-axis limits
  p <- p + coord_cartesian(xlim = xLim, ylim = yLim)
  
  # if (!is.null(margin)){p <- p + theme(plot.margin = unit(margin, "cm"))}
  
  # -------------------------------------------------------------------------- #
  ### Close saving:
  
  if (isSave){
    print(p)
    dev.off()
  }
  print(p)
  return(p)
  cat("Finished :-)\n")
}

# ============================================================================ #
#### Clean formula for new lines and excessive whitespaces: ####

clean_formula <- function(formula){
  #' Clean formula split over multiple rows by removing newline characters and double whitespaces.
  #' @param formula scalar string, formula to clean.
  #' @return output scalar string, cleaned-up formula string.
  
  output <- formula
  output <- gsub("\n", "", output)
  output <- gsub("  ", "", output)
  
  return(output)
  
}

# ============================================================================ #
#### Fit generalized additive mixed model (GAMM) based on input variables: ####

fit_gamm <- function(data, yVar, timeVar, splitVar, groupVar, bsType, addARIMA = NULL, useScat = F, isPrint = T){
  #' Retrieve previously fitted or fit anew a generalized additive mixed-effects model (GAMM)
  #' based on input variables.
  #' @param data data frame, trial-by-trial data.
  #' @param yVar scalar string, name of dependent variable (must be numeric).
  #' @param timeVar scalar string, name of time-course variable (must be numeric).
  #' @param splitVar scalar string, name of condition variable to split by (must be factor).
  #' @param groupVar scalar string, name of variable by which data is grouped (must be factor).
  #' @param bsType scalar string, type of random effects, either
  #' - "fe" = only fixed effects, no random effects.
  #' - "ri" = only random intercept.s
  #' - "rs" = random intercepts and slopes.
  #' - "fs" = random smooths.
  #' @param addARIMA Boolean, use AR(1) model (default: F for fe/ri/fs, T, for fs). 
  #' @param useScat Boolean, use scaled t model (default: F).
  #' @return mod fitted model object. 
  
  ## Print specifics to console:
  cat(paste0("Fit GAMM with settings DV = ", yVar, ", time = ", timeVar, ", split by = ", splitVar, ", grouped by ", groupVar, "\n"))
  cat(paste0("Random effects are of type ", bsType, "\n"))
  
  # ------------------------------------------------------------------------ #
  ## Check inputs:
  
  if(!is.numeric(data[, yVar])){stop("yVar must be numeric")}
  if(!is.numeric(data[, timeVar])){stop("timeVar must be numeric")}
  if(!is.factor(data[, splitVar])){stop("splitVar must be a factor")}
  if(!is.factor(data[, groupVar])){stop("groupVar must be a factor")}
  
  ## Use AR(1) model as default for FS:
  if(is.null(addARIMA)){if(bsType == "fs"){addARIMA <- T} else {addARIMA <- F}}
  if(addARIMA){cat("Add AR(1) model\n")}
  if(useScat){cat("Use scaled t model for heavy tailed response variables\n")}
  
  # ------------------------------------------------------------------------ #
  ### Determine model name:
  
  modName <- paste0("gamm_", yVar, "~", timeVar, "_x_", splitVar, "_per_", groupVar, "_", bsType)
  if(addARIMA){  modName <- paste0(modName, "_AR")}
  if(useScat){  modName <- paste0(modName, "_scat")}
  modName <- paste0(modName, ".rds")
  
  # ------------------------------------------------------------------------ #
  ### Check if already exists; if yes, retrieve; if not, fit anew:
  
  if (file.exists(paste0(dirs$modelDir, modName))){
    
    cat(paste0(">>> Found ", modName, ", load \n"))
    mod <- readRDS(paste0(dirs$modelDir, modName))
    
  } else {
    
    # ---------------------------------------------------------------------- #
    ### Specify formula:
    
    if (bsType  == "fe"){ # if only fixed effects
      
      # do not include main effect of time
      formula <- paste0(yVar, " ~ ", splitVar, " + 
                        s(", timeVar, ") + 
                        s(", timeVar, ", by = ", splitVar, ")")
      
    } else if (bsType == "ri"){ # if random intercepts
      
      # do not include main effect of time
      formula <- paste0(yVar, " ~ ", splitVar, " + 
                  s(", timeVar, ", by = ", splitVar, ") + 
                  s(", groupVar, ", bs = \"re\")")
      
    } else if (bsType == "rs"){ # if random intercepts and random slopes
      
      # do not include main effect of time
      formula <- paste0(yVar, " ~ ", splitVar, " + 
                  s(", timeVar, ", by = ", splitVar, ") + 
                  s(", groupVar, ", bs = \"re\") + 
                  s(", groupVar, ", ", timeVar, ", bs = \"re\")")
      
    } else if (bsType == "fs"){ # if random smooths
      
      # do not include main effect of time
      # random smooths covers for both intercept and slope (it's just a smooth, no separation into intercept and slope)
      formula <- paste0(yVar, " ~ ", splitVar, " + 
                  s(", timeVar, ", by = ", splitVar, ") + 
                  s(", timeVar, ", ", groupVar, ", bs = \"fs\", m = 1)")
      
    } else { # if unknown
      
      stop(paste0("bsType = ", bsType, " is unknown"))
      
    } # end of formula loop
    
    # ---------------------------------------------------------------------- #
    ## Add lag-1 correlation for AR(1) model:
    if(addARIMA){
      
      ## Add start of new time series to data:
      data <- start_event(data, column = timeVar, event = groupVar)
      
      ## Estimate auto-correlation from random-slopes model above as starting value:
      if(exists("mod")){
        (valRho <- acf(resid(mod), plot = FALSE)$acf[2]) # estimate auto-correlation 
        cat(paste0("Use auto-correlation value of valRho = ", valRho, " estimated from last model\n"))  
      } else {
        valRho <- -0.03619877
        cat(paste0("No past mod to estimate auto-correlation; set valRho = ", valRho, "\n"))
      }
      
    }
    
    # ---------------------------------------------------------------------- #
    ### Fit model:
    
    ## Start time:
    start.time <- Sys.time();
    cat(paste0(">>> Fit generalized additive mixed-effects model (GAMM) with formula\n", clean_formula(formula), "\n"))
    cat("... \n")
    
    if(addARIMA & useScat){
      mod <- bam(formula = eval(parse(text = formula)),
                 AR.start = data$start.event, rho = valRho,
                 family = "scat",
                 data = data, discrete = T) #
    } else if(addARIMA){
      mod <- bam(formula = eval(parse(text = formula)),
                 AR.start = data$start.event, rho = valRho,
                 data = data, discrete = T) #
    } else if(useScat){
      mod <- bam(formula = eval(parse(text = formula)),
                 family = "scat",
                 data = data, discrete = T) #
    } else {
      mod <- bam(formula = eval(parse(text = formula)),
                 data = data, discrete = T) # 
    }
    
    ## Stop time:
    end.time <- Sys.time(); beep()
    dif <- difftime(end.time,start.time); print(dif)
    
    # ------------------------------------------------------------------------ #
    ### Save model:
    
    cat(paste0(">>> Save ", modName, "\n"))
    saveRDS(mod, paste0(dirs$modelDir, modName))
    cat(">>> Saved :-)\n")
    
  }
  
  # -------------------------------------------------------------------------- #
  ### Print output:
  
  if(isPrint){
    print(summary(mod))
    # report_stats(mod)
  }
  cat("Finished! :-)\n")
  
  return(mod)
  
}

# ============================================================================ #
#### Fit generalized additive mixed model (GAMM) to test condition difference: ####

fit_gamm_diff <- function(data, yVar, timeVar, splitVar, splitLevels = NULL, groupVar, bsType, addARIMA = NULL, useScat = F, isPrint = F){
  #' Retrieve previously fitted or fit anew a generalized additive mixed-effects model (GAMM)
  #' based on input variables.
  #' @param data data frame, trial-by-trial data.
  #' @param yVar scalar string, name of dependent variable (must be numeric).
  #' @param timeVar scalar string, name of time-course variable (must be numeric).
  #' @param splitVar scalar string, name of condition variable to split by (must be factor).
  #' @param splitLevels vector of 2 strings, levels of splitVar to retain (optional.
  #' @param groupVar scalar string, name of variable by which data is grouped (must be factor).
  #' @param bsType scalar string, type of random effects, either
  #' - "fe" = only fixed effects, no random effects.
  #' - "ri" = only random intercept.s
  #' - "rs" = random intercepts and slopes.
  #' - "fs" = random smooths.
  #' @param addARIMA Boolean, use AR(1) model (default: F for fe/ri/fs, T, for fs). 
  #' @param useScat Boolean, use scaled t model (default: F).
  #' @return mod fitted model object. 
  
  ## Print specifics to console:
  cat(paste0("Fit GAMM with settings DV = ", yVar, ", time = ", timeVar, ", split by = ", splitVar, ", grouped by ", groupVar, "\n"))
  cat(paste0("Random effects are of type ", bsType, "\n"))
  
  # ------------------------------------------------------------------------ #
  ## Check inputs:
  
  if(!is.numeric(data[, yVar])){stop("yVar must be numeric")}
  if(!is.numeric(data[, timeVar])){stop("timeVar must be numeric")}
  if(!is.factor(data[, splitVar])){stop("splitVar must be a factor")}
  if(!is.factor(data[, groupVar])){stop("groupVar must be a factor")}
  
  ## Use AR(1) model as default for FS:
  if(is.null(addARIMA)){if(bsType == "fs"){addARIMA <- T} else {addARIMA <- F}}
  if(addARIMA){cat("Add AR(1) model\n")}
  if(useScat){cat("Use scaled t model for heavy tailed response variables\n")}
  
  # ------------------------------------------------------------------------ #
  ### Retrieve and check levels to splitVar:
  
  if (is.null(splitLevels))
    splitLevels <- levels(data[, selVar])
  end
  if(!length(splitLevels)){stop("splitLevels must have 2 elements")}
  
  # ------------------------------------------------------------------------ #
  ### Select data:
  
  cat(paste0("Select levels ", paste0(splitLevels, collapse = ", "), " of variable ", splitVar, "\n"))
  selData <- droplevels(data[which(data[, splitVar] %in% splitLevels), ])
  cat(paste0("Test difference ", paste0(splitLevels, collapse = " - "), "\n"))
  
  # ------------------------------------------------------------------------ #
  ### Create difference variable:
  
  # splitVarOF <- paste0(splitVar, "_OF")
  # selData[, splitVarOF] <- as.factor(selData[, splitVar]) # factor
  # selData[, splitVarOF] <- as.ordered(selData[, splitVarOF]) # ordered factor
  # contrasts(selData[, splitVarOF]) <- "contr.treatment"
  
  splitVarOF <- "difference_OF"
  selData$difference_OF <- as.factor(selData[, splitVar]) # factor
  selData$difference_OF <- as.ordered(selData$difference_OF) # ordered factor
  contrasts(selData$difference_OF) <- "contr.treatment"
  
  # ------------------------------------------------------------------------ #
  ### Determine model name:
  
  modName <- paste0("gamm_diff_", yVar, "~", timeVar, "_x_", splitVar, "_",
                    paste0(splitLevels, collapse = "_"),
                    "_per_", groupVar, "_", bsType)
  if(addARIMA){modName <- paste0(modName, "_AR")}
  if(useScat){modName <- paste0(modName, "_scat")}
  modName <- paste0(modName, ".rds")
  
  # ------------------------------------------------------------------------ #
  ### Check if already exists; if yes, retrieve; if not, fit anew:
  
  if (file.exists(paste0(dirs$modelDir, modName))){
    
    cat(paste0(">>> Found ", modName, ", load \n"))
    mod <- readRDS(paste0(dirs$modelDir, modName))
    
  } else {
    
    # ---------------------------------------------------------------------- #
    ### Specify formula:
    
    if (bsType  == "fe"){ # if only fixed effects
      
      # do not include main effect of time
      formula <- paste0(yVar, " ~ ", splitVarOF, " + 
                  s(", timeVar, ") + 
                  s(", timeVar, ", by = ", splitVarOF, ")")
      
    } else if (bsType == "ri"){ # if random intercepts
      
      formula <- paste0(yVar, " ~ ", splitVarOF, " + 
                  s(", timeVar, ") + 
                  s(", timeVar, ", by = ", splitVarOF, ") + 
                  s(", groupVar, ", bs = \"re\")")
      
    } else if (bsType == "rs"){ # if random intercepts and random slopes
      
      formula <- paste0(yVar, " ~ ", splitVarOF, " + 
                  s(", timeVar, ") + 
                  s(", timeVar, ", by = ", splitVarOF, ") + 
                  s(", groupVar, ", bs = \"re\") + 
                  s(", groupVar, ", ", timeVar, ", bs = \"re\")")
      
      
    } else if (bsType == "fs"){ # if random smooths
      
      # random smooths covers for both intercept and slope (it's just a smooth, no separation into intercept and slope)
      formula <- paste0(yVar, " ~ ", splitVarOF, " + 
                  s(", timeVar, ") + 
                  s(", timeVar, ", by = ", splitVarOF, ") + 
                  s(", timeVar, ", ", groupVar, ", bs = \"fs\", m = 1)")
      
    } else { # if unknown
      
      stop(paste0("bsType = ", bsType, " is unknown"))
      
    } # end of formula loop
    
    
    # ---------------------------------------------------------------------- #
    ## Add lag-1 correlation for AR(1) model:
    
    if(addARIMA){
      
      ## Add start of new time series to data:
      data <- start_event(data, column = timeVar, event = groupVar)
      
      ## Estimate auto-correlation from random-slopes model above as starting value:
      if(exists("mod")){
        (valRho <- acf(resid(mod), plot = FALSE)$acf[2]) # estimate auto-correlation 
        cat(paste0("Use auto-correlation value of valRho = ", valRho, " estimated from last model\n"))  
      } else {
        valRho <- -0.03619877
        cat(paste0("No past mod to estimate auto-correlation; set valRho = ", valRho, "\n"))
      }
    }
    
    # ---------------------------------------------------------------------- #
    ### Fit model:
    
    ## Start time:
    start.time <- Sys.time();
    cat(paste0(">>> Fit generalized additive mixed-effects model (GAMM) with formula\n", clean_formula(formula), "\n"))
    cat("... \n")
    
    if(addARIMA & useScat){
      mod <- bam(formula = eval(parse(text = formula)),
                 AR.start = selData$start.event, rho = valRho,
                 family = "scat",
                 data = selData, discrete = T) #
    } else if(addARIMA){
      mod <- bam(formula = eval(parse(text = formula)),
                 AR.start = selData$start.event, rho = valRho,
                 data = selData, discrete = T) #
    } else if(useScat){
      mod <- bam(formula = eval(parse(text = formula)),
                 family = "scat",
                 data = selData, discrete = T) #
    } else {
      mod <- bam(formula = eval(parse(text = formula)),
                 data = selData, discrete = T) # 
    }    
    ## Stop time:
    end.time <- Sys.time(); beep()
    dif <- difftime(end.time,start.time); print(dif)
    
    # ------------------------------------------------------------------------ #
    ### Save model:
    
    cat(paste0(">>> Save ", modName, "\n"))
    saveRDS(mod, paste0(dirs$modelDir, modName))
    cat(">>> Saved :-)\n")
    
  }
  
  # -------------------------------------------------------------------------- #
  ### Print output:
  
  if (isPrint){
    print(summary(mod))
    # report_stats(mod)
  }
  cat("Finished! :-)\n")
  
  return(mod)
  
}

# ============================================================================ #
#### Compute CIs given mgcv model outut: ####

quickCImgcv <- function(mod, level = 0.95, nRound = 3){
  #' Compute CIs for given GAMM model fit with mgcv given SEs from model.
  #' @param data mod model objected fitted with lme4.
  #' @param level numeric, 0-1, level of CIs (default: 0.95).
  #' @param nRound integer, number of digits after comma to round to (default: 2).
  #' @return print to console.
  
  ## Checks:
  if((class(mod)[1] != "gam") & (class(mod)[1] != "bam")){stop("Provided model is not of class gam")}
  
  ## Determine z-value by which to multiply SE:
  twoSideLevel <- 1 - (1 - level) / 2 # correct for two-sided test
  zVal <- qnorm(twoSideLevel) # respective z-level threshold
  
  ## Predictor index:
  # regIdx <- which(grepl("difference_OF", names(coef(mod)), fixed = T))[1] # index of difference variable
  regIdx <- 2
  
  cat(paste0("Compute " , level*100, "%-CIs for predictor \"", names(coef(mod))[regIdx], "\"... \n"))
  
  # -------------------------------------------------------------------------- #
  ### Compute CIs:
  # https://stats.stackexchange.com/questions/364568/extracting-significance-of-gam-smoothing-terms
  tmp <- summary.gam(mod) # extract summary
  # names(mod$model)
  
  ## Regression coefficient:
  # fitB <- as.numeric(coef(mod)[regIdx]) # estimate
  fitB <- as.numeric(tmp$p.coeff[regIdx])
  if (fitB < 0){fitB <- fitB * -1} # invert if negative
  
  ## Standard error:
  # fitSE <- as.numeric(sqrt(diag(vcov(mod, unconditional  = TRUE)))[regIdx]) # se
  fitSE <- as.numeric(tmp$se[regIdx])
  
  ## Print to console:
  cat(paste0("b = ", round(fitB, nRound), ", 95%-CI [", 
             round(fitB - zVal * fitSE, nRound), ", ",
             round(fitB + zVal * fitSE, nRound), "], "))
  
  # -------------------------------------------------------------------------- #
  ## Print t-value:
  # https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam
  # See mgc manual https://cran.r-project.org/web/packages/mgcv/mgcv.pdf
  # documentation for function summary.gam, page 295.
  cat(paste0("t(", round(tmp$residual.df), ") = ", round(tmp$p.t[regIdx], nRound), 
             ", p = ", round(tmp$p.pv[regIdx], nRound), "\n"))
  
}

# ============================================================================ #
#### Plot all smooths from mgcv model with thicker lines: ####

custom_plot_smooth <- function(mod,
                               colVec = NULL, ltyVec = NULL, SEweight = 1, addLegend = F,
                               xLab = NULL, yLab = NULL, main = NULL, isPNG = T){
  
  # ------------------------------------------------------------------------ #
  ### Fixed settings:
  
  if (isPNG){
    CEX <- 2.5
    LWD <- 2
    labelInt <- 3.5
    labelSlope <- 0.65
    xMGP <- c(2.5, 1.5, 0)
    yMGP <- c(2, 1, 0)
  } else { # without saving:
    CEX <- 1.5
    LWD <- 3
    labelInt <- 2.5
    labelSlope <- 0.5
    xMGP <- c(1.5, 1, 0)
    yMGP <- c(2, 1, 0)
  }  
  
  # ------------------------------------------------------------------------ #
  ### Recover variable names:
  
  varNames <- names(mod$model)
  yVar <- varNames[1]
  splitVar <- varNames[2]
  timeVar <- varNames[3]
  groupVar <- varNames[4]
  
  # ------------------------------------------------------------------------ #
  ### Retrieve plotting settings:
  
  xLim <- c(0, max(mod$model[, timeVar]))
  if(is.null(colVec)){colVec <- retrieve_colour(splitVar)}
  if(is.null(ltyVec)){
    nLevel <- length(levels(data[, splitVar]))
    ltyVec <- rep(1, nLevel)
  }
  if(is.null(xLab)){xLab <- substitute_label(timeVar)}
  if(is.null(yLab)){yLab <- substitute_label(yVar)}
  # if(is.null(main)){main <- ""}
  
  if(addLegend){
    yMean <- as.numeric(tapply(data[, yVar], interaction(data[, timeVar], data[, splitVar]), mean, na.rm = T))
    legendPos <- c(xLim[2] * 0.90, min(yMean) + (max(yMean) - min(yMean))*.83)
    cat(paste0("Plot legend at position ", legendPos[1], ", ", legendPos[2], "\n"))
  } else {
    legendPos <- c(xLim[2] + 1, 0)
  }
  
  # ------------------------------------------------------------------------ #
  ### Save:
  
  if(isPNG){
    
    if(length(dev.list()!=0)){dev.off()} # close any open plots
    
    if(any(grepl("bs = \"fs\"", mod$formula))){ 
      bsType <- "fs"
    } else if (any(grepl("bs = \"re\"", mod$formula))){
      bsType <- "rs"
    } else {
      bsType <- "fe"
    }
    
    plotName <- paste0("gam_", yVar, "~", timeVar, "_x_", splitVar, "_",
                       "_per_", groupVar, "_", bsType)
    if(mod$AR1.rho != 0){plotName <- paste0(plotName, "_AR")}
    if(grepl("Scaled t", mod$family$family)){plotName <- paste0(plotName, "_scat")}
    if(addLegend){plotName <- paste0(plotName, "_addLegend")}
    cat(paste0("Save plot under ", plotName, ".png\n"))
    png(paste0(dirs$plotDir, plotName, ".png"), width = 600, height = 480)
    
  }
  
  # ------------------------------------------------------------------------ #
  ### Plot all levels:
  
  ## Change margin to have more space on the left:
  if (isPNG){
    par(mar = c(5.1, 7.5, 1.1, 0)) # c(bottom, left, top, right))
  } else {
    par(mar = c(5.1, 5.6, 1.1, 0)) # c(bottom, left, top, right))
  }
  
  p <- plot_smooth(mod, view = timeVar, plot_all = splitVar, 
                   col = colVec, lty = ltyVec, lwd = 6, # 3
                   cex.main = CEX, cex.axis = CEX, cex.lab = CEX,
                   xlim = xLim, las = 1, h0 = NULL, 
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main,
                   # legend_plot_all = "topright",
                   legend_plot_all = list(x = legendPos[1], y = legendPos[2]),
                   # main = paste0("All levels of ", tolower(substitute_label(splitVar))), 
                   se = SEweight, rm.ranef = TRUE, rug = FALSE)
  
  # ------------------------------------------------------------------------ #
  ### Make axes and highlight lines thicker:
  
  # p <- plot_smooth(mod, view = timeVar, plot_all = splitVar)
  
  ## Retrieve y-axis labels:
  # yLim <- c(ceiling(min(p$fv$fit - p$fv$CI)), floor(max(p$fv$fit + p$fv$CI)))
  # yLim <- c(ceiling(min(p$fv$fit - p$fv$CI)), floor(max(p$fv$fit + p$fv$CI)))
  # yLim <- c(round(min(p$fv$fit - p$fv$CI)), round(max(p$fv$fit + p$fv$CI)))
  yLim <- c(min(p$fv$fit - p$fv$CI), max(p$fv$fit + p$fv$CI))
  yLim <- find_round_lim(yLim)
  
  ## Thicker axes:
  xStep <- find_step(xLim, nTickTarget = 5) # tick step size x-axis
  yStep <- find_step(yLim, nTickTarget = 5) # tick step size y-axis
  yDigit <- 1 + abs(floor(log10(yStep))) + as.numeric(abs(yStep) < 1) # number digits y-axis ticks
  # yMGP[1] <- yMGP[1] + yDigit - 1
  cat(paste0("Use x-axis limits ", xLim[1], "-", xLim[2], " in steps of ", xStep, "\n"))
  axis(1, at = seq(xLim[1], xLim[2], xStep), lwd = LWD, lwd.tick = LWD, cex.axis = CEX, lab = T, mgp = xMGP)
  cat(paste0("Use y-axis limits ", yLim[1], "-", yLim[2], " in steps of ", yStep, "\n"))
  axis(2, at = seq(yLim[1], yLim[2], yStep), lwd = LWD, lwd.tick = LWD, cex.axis = CEX, las = 1, lab = T, mgp = yMGP)
  title(xlab = xLab, line = labelInt, cex.lab = CEX)
  title(ylab = yLab, line = labelInt + (yDigit - 1) * labelSlope, cex.lab = CEX)
  
  if(isPNG){dev.off()}
  
  ## Change margin back to default:
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))
  
  cat("Finished! :-)")
  
}

# ============================================================================ #
#### Plot difference term from mgcv model with thicker lines: ####

custom_plot_diff <- function(mod, 
                             selCol = "black", SEweight = 1.96,
                             xLab = NULL, yLab = NULL, main = NULL, isPNG = T){
  
  # ------------------------------------------------------------------------ #
  ### Fixed settings:
  
  if (isPNG){
    CEX <- 2.5
    LWD <- 2
    labelInt <- 3.5
    labelSlope <- 0.65
    xMGP <- c(2.5, 1.5, 0)
    yMGP <- c(2, 1, 0)
  } else { # without saving:
    CEX <- 1.5
    LWD <- 3
    labelInt <- 2.5
    labelSlope <- 0.5
    xMGP <- c(1.5, 1, 0)
    yMGP <- c(2, 1, 0)
  }  
  
  # ------------------------------------------------------------------------ #
  ### Recover variable names:
  
  varNames <- names(mod$model)
  yVar <- varNames[1]
  if(!exists("splitVar")){splitVar <- varNames[2]} # use global variable if possible
  timeVar <- varNames[3]
  groupVar <- varNames[4]
  
  # ------------------------------------------------------------------------ #
  ### Retrieve plotting settings:
  
  xLim <- c(0, max(mod$model[, timeVar]))
  
  if(is.null(xLab)){xLab <- substitute_label(timeVar)}
  if(is.null(yLab)){yLab <- paste0("Est. diff. in ", substitute_label(yVar))} # tolower()
  if(is.null(main)){main <- ""}
  
  # ------------------------------------------------------------------------ #
  ### Save:
  
  if(isPNG){
    
    if(length(dev.list()!=0)){dev.off()} # close any open plots
    
    if(any(grepl("bs = \"fs\"", mod$formula))){ 
      bsType <- "fs"
    } else if (any(grepl("bs = \"re\"", mod$formula))){
      bsType <- "rs"
    } else {
      bsType <- "fe"
    }
    
    plotName <- paste0("gam_diff_", yVar, "~", timeVar, "_x_", splitVar, "_",
                       paste0(splitLevels, collapse = "_"),
                       "_per_", groupVar, "_", bsType)
    if(mod$AR1.rho != 0){plotName <- paste0(plotName, "_AR")}
    if(grepl("Scaled t", mod$family$family)){plotName <- paste0(plotName, "_scat")}
    cat(paste0("Save plot under ", plotName, ".png\n"))
    png(paste0(dirs$plotDir, plotName, ".png"), width = 600, height = 480)
    
  }
  
  # ------------------------------------------------------------------------ #
  ### Plot difference:
  
  ## Change margin to have more space on the left:
  par(mar = c(5.1, 5.6, 1.1, 0)) # c(bottom, left, top, right))
  
  # p <- plot_diff(mod, view = timeVar, comp = list(difference_OF = c(splitLevels[1], splitLevels[2])))
  p <- plot_diff(mod, view = timeVar, comp = list(difference_OF = c(splitLevels[1], splitLevels[2])), 
                 xlim = xLim, # transform.view = NULL,
                 col = selCol, lwd = 3, las = 1, # h0 = NULL, 
                 cex.lab = CEX, cex.axis = CEX, cex.main = CEX,
                 col.diff = plotfunctions::alpha("red", f = 0), # switch off
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main,
                 # main = paste0(splitLevels[1], " minus ", splitLevels[2]), 
                 rm.ranef = F, add = F) # 
  par(ask = F)
  
  # ------------------------------------------------------------------------ #
  ### Make axes and highlight lines thicker:
  
  ## Retrieve y-axis labels:
  # p <- plot_diff(mod, view = timeVar, comp = list(difference_OF = c(splitLevels[1], splitLevels[2])))
  # yLim <- c(ceiling(min(p$est - p$CI)), floor(max(p$est - p$CI)))
  # yLim <- c(ceiling(min(p$est - p$CI)), floor(max(p$est - p$CI)))
  # yLim <- c(round(min(p$est - p$CI)), round(max(p$est - p$CI)))
  yLim <- c(min(p$est - p$CI), max(p$est + p$CI))
  yLim <- find_round_lim(yLim)
  
  ## Thicker axes:
  xStep <- find_step(xLim, nTickTarget = 5) # tick step size x-axis
  yStep <- find_step(yLim, nTickTarget = 5) # tick step size y-axis
  yDigit <- 1 + abs(floor(log10(yStep))) + as.numeric(abs(yStep) < 1) # number digits y-axis ticks
  # yMGP[1] <- yMGP[1] + yDigit - 1
  axis(1, at = seq(xLim[1], xLim[2], xStep), lwd = LWD, lwd.tick = LWD, cex.axis = CEX, lab = T, mgp = xMGP)
  axis(2, at = seq(yLim[1], yLim[2], yStep), lwd = LWD, lwd.tick = LWD, cex.axis = CEX, las = 1, lab = T, mgp = yMGP)
  title(xlab = xLab, line = labelInt, cex.lab = CEX)
  title(ylab = yLab, line = labelInt + (yDigit - 1) * labelSlope, cex.lab = CEX)
  
  
  ## Thicker highlight lines:
  x <- find_difference(p$est, p$CI, f = 1, xVals = p[, timeVar])
  if(!is.null(x)){
    xStart <- x$start
    xStop <- x$end
    for(iElem in 1:length(xStart)){
      cat(paste0("Significant difference from ", xStart[iElem], " - ", xStop[iElem], "\n"))
      axis(1, at = c(xStart[iElem], xStop[iElem]), lwd = LWD*2, col = "red", labels = F)
      abline(v = xStart[iElem], col = "red", lwd = LWD, lty = 2)
      abline(v = xStop[iElem], col = "red", lwd = LWD, lty = 2)
    }
  }
  
  if(isPNG){dev.off()}
  
  ## Change margin back to default:
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))
  
  cat("Finished! :-)")
  
}

# END
