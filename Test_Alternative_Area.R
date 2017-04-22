#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Test Alternative Areas/Regions
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.28.17
#
#Purpose: Example of How to Use helper functions
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
source('R/create-VAST-input.r')
source('R/plot-VAST-output.r')

require(VAST)
require(TMB)

#=======================================================================
##### SETUP INPUT DATA #####

#Generate a dataset
species.codes<- 21740 #c(30420) #Rockfish
lat_lon.def="mean"