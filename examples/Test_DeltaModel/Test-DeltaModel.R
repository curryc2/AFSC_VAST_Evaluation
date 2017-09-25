#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare VAST estimates with delta-GLMM
#
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 9.13.17
#
#Purpose: Evaluate whether differences between design-based and VAST (model-based) indices and uncertainty, is a function of the delta model structure. 
#
#
#==================================================================================================
#NOTES:
#  a) Spiny Dogfish removed from evaluation because design based index is 0 for Western GOA in 2013
#  b) Non-convergence with intercepts estimted as IaY i.e. ==1, independent of spatio-temporal RE specs.
#==================================================================================================
#TIMING:

##==================================================================================================

require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)
require(R2admb)
require(reshape2)
require(gridExtra)
require(ggthemes)
require(cowplot)