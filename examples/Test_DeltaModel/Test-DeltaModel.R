#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Apportionment Between VAST and
#                                                              RE Model for GOA
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 5.11.17
#
#Purpose: Evaluate apportionment (Eastern, Central, Western) for Gulf of Alaska RACE
#           Bottom Trawl indices of abundance. 
#             1) Evaluate across three alternative knot numbers
#             2) Several specifications for autocorrelation in EC and PCR intercept
#
#
#==================================================================================================
#NOTES:
#  a) Spiny Dogfish removed from evaluation because design based index is 0 for Western GOA in 2013
#  b) Non-convergence with intercepts estimted as IaY i.e. ==1, independent of spatio-temporal RE specs.
#==================================================================================================
#TIMING:
# [1] "### START: Tue Jul 18 15:32:25 2017"
# [1] "### END: Wed Jul 19 01:35:23 2017"
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