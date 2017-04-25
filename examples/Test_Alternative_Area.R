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
source('R/calc-design-based-index.r')
source('R/create-Data-Geostat.r')
source('R/load-RACE-data.r')

require(VAST)
require(TMB)
require(ggmap)
require(ggplot2)

#=======================================================================
##### SETUP INPUT DATA #####

#Generate a dataset
species.codes <- 21720 #Pacific Cod #c(30420) #Rockfish

#DERIVED OBJECTS
Region  <- "Gulf_of_Alaska"


dat.goa <- load_RACE_data(species.codes=species.codes, area='GOA')

dat.ai <- load_RACE_data(species.codes=species.codes, area='AI')

dat.bs <- load_RACE_data(species.codes=species.codes, area='BS')

#Test next level function
dgs.goa <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="GOA")
dgs.ai <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="AI")

#
ggmap()

#Lets get a map
map.input <- dat.bs

#Left, bottom, right, top
map.dat <- get_map(location=c(min(map.input$START_LONGITUDE),min(map.input$START_LATITUDE),max(map.input$START_LONGITUDE), max(map.input$START_LATITUDE)),
                     maptype='terrain')
ggmap(map.dat)
