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
#
#  
#==================================================================================================
#TIMING:
#
#==================================================================================================

require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Autoregressive")


#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

#Limit species included
species.list <- species.list[species.list$include=='Y',]
n.species <- nrow(species.list)

#Create
species.series <- c(1:n.species)

#=======================================================================
##### CONTROL SECTION #####
#Number of cores to use
n.cores <- detectCores()-1

#Boolean for running estimation models
do.estim <- TRUE

#Trial Knot Numbers
trial.knots <- c(100,500,1000)
n.trial.knots <- length(trial.knots)

#Trial AUTOREGRESSIVE specifications
#Note starts at 0
# rho.int.types <- c('Fixed_Effect','Independent_Among_Years','Random_Walk','Constant_Intercept','Autoregressive')
rho.int.types <- c('FE','IaY','RW','CI','AR')
# rho.stRE.types <- c('Independent_Among_Years',NA,'Random_Walk',NA,'Autoregressive')
rho.stRE.types <- c('IaY',NA,'RW',NA,'AR')

#Read in Autoregressive Input
trial.rho <- t(read.csv('Data/Test-Autoregressive-Input.csv', header=TRUE, stringsAsFactors=FALSE)[,-c(1:2)])
n.trial.rho <- nrow(trial.rho)

# #Intercept
# rho.inter.ep <- c(0) #Encounter Probability Component 
# rho.inter.pcr <- c(0) #Positive Catch Rate Component
# #Spatio-temporal RE
# rho.stRE.ep #Encounter Probability Component 
# rho.stRE.pcr #Positive Catch Rate Component

#Boolean for bias correction
bias.correct <- FALSE
#=======================================================================
##### Run VAST model  #####
Version <- "VAST_v2_4_0"
lat_lon.def <- "start"

#SPATIAL SETTINGS
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
# n_x = c(100, 250, 500, 1000, 2000)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"))