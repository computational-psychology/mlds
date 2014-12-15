library(MLDS)

DefineMyScale <- function(rr = c(seq(0, 0.9, len = 10), 0.98)) {
  #sqrt(r)  # for equal-spacing in r^2
  rr
}


stim <- do.call(DefineMyScale , list())

NumStimuli <- length(stim)

allTrials <- t(combn(seq(NumStimuli), 4))  # 3 for triads, 4 for quadruples

# randomize order of trials
trialOrder <- sample(seq(nrow(allTrials)), replace = FALSE)
trial <- allTrials[trialOrder, ]

# randomize pairs left or right 
topbot <- as.logical(rbinom(nrow(trial), 1, 0.5))
trial[topbot, ] <- trial[topbot, c(3, 4, 1, 2)]    # c(3, 2, 1)   for triads

# empty vector for responses
resp <- rep(NA, nrow(trial))


