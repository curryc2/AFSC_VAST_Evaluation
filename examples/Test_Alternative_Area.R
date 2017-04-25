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

#==============================================================================
#TEST: reading RACE data
dat.goa <- load_RACE_data(species.codes=species.codes, area='GOA')

dat.ai <- load_RACE_data(species.codes=species.codes, area='AI')

dat.ebs <- load_RACE_data(species.codes=species.codes, area='EBS') #DOESN'T WORK "bs" is not an area in RACE data
dat.bs <- load_RACE_data(species.codes=species.codes, area='BS') #Note that race data does not separate EBS and BS

#==============================================================================
#TEST: creating Data_Geostat object
dgs.goa <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="GOA")
dgs.ai <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="AI")

dgs.ebs <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="EBS") #Doesn't work!!!

dgs.bs <- create_Data_Geostat(species.codes, lat_lon.def="mean", area="BS")

#==============================================================================
#TEST: creating VAST input
# This is the only function that takes Region in place of area.
# It gets converted right off the bat.
# 
# create_Data_Geostat() uses "area"
# Prepare_Extrapolation_Data_Fn() uses Region

#NOTE: Might consider changing the input type around cor create VAST input

#==============================================================================
#TEST: Design-based Index Function




#
ggmap()

#Lets get a map
map.input <- dgs.goa



#Left, bottom, right, top
# map.dat <- get_map(location=c(lon = mean(dgs.goa$Lon), lat = mean(dgs.goa$Lat)),
#                      maptype='terrain', source='stamen', crop=FALSE)

#Left, bottom, right, top
map.dat <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                              max(map.input$Lon),max(map.input$Lat)),
                   maptype='terrain', source='stamen', crop=FALSE)

map.dat.1 <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                                max(map.input$Lon),max(map.input$Lat)),
                     maptype='watercolor', source='stamen', crop=FALSE)

map.dat.2 <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                                max(map.input$Lon),max(map.input$Lat)),
                     maptype='toner-lite', source='stamen', crop=FALSE)

map.dat.3 <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                                max(map.input$Lon),max(map.input$Lat)),
                     maptype='toner', source='stamen', crop=FALSE)

#DON'T USE
# map.dat.4 <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
#                                 max(map.input$Lon),max(map.input$Lat)),
#                      maptype='terrain', source='google', crop=FALSE)

ggmap(map.dat) #terrain: staimen
# ggmap(map.dat.4) #terrain: google
ggmap(map.dat.1) #watercolor
ggmap(map.dat.2) #toner-lite
ggmap(map.dat.3) #toner





g <- ggmap(map.dat.2) +
       geom_point(data=map.input, aes(x=Lon, y=Lat, color=Year), alpha=1, size=0.25)

g


#Plot bs
map.input <- dgs.bs

map.dat <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                                max(map.input$Lon),max(map.input$Lat)),
                                maptype='toner-lite', source='stamen', crop=FALSE)

g.bs <- ggmap(map.dat) +
          geom_point(data=map.input, aes(x=Lon, y=Lat, color=Year),
                       alpha=1, size=0.25)
g.bs

#Plot goa
map.input <- dgs.goa

map.dat <- get_map(location=c(min(map.input$Lon),min(map.input$Lat),
                              max(map.input$Lon),max(map.input$Lat)),
                   maptype='toner-lite', source='stamen', crop=FALSE)

g.goa <- ggmap(map.dat) +
          geom_point(data=map.input, aes(x=Lon, y=Lat, color=Year),
                       alpha=1, size=0.25)
g.goa

#Plot AI
#Plot bs
# map.input <- dgs.ai
# 
# map.dat <- get_map(location=c(179,min(map.input$Lat),
#                               -150,max(map.input$Lat)),
#                    maptype='toner-lite', source='stamen', crop=FALSE)
# 
# g.ai <- ggmap(map.dat) +
#   geom_point(data=map.input, aes(x=Lon, y=Lat, color=Year),
#              alpha=1, size=0.25)
# g.ai