#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Testing Observation Model Specification
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.30.17
#
#Purpose: To test impact of changing distribution on observation model for positive catch rates,
#           across species.
#  1) Loop through obs model specifications for positive catch rate distribution.
#  2) 
#
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
require(snowfall)
require(parallel)
require(ggplot2)
require(TMB)
require(TMBhelper)
require(VAST)

#Source necessary files
source("R/create-VAST-input.r")

trial.obs <- c(0:7)
n.trial.obs <- length(trial.obs)

#Names for trial Obs models
names.trial.obs <- c('Normal','Lognormal','Gamma','Negative_binomial','Conway-Maxwell_Poisson','Poisson')




t <- 1
for(t in 1:n.trial.obs) {
  
}
