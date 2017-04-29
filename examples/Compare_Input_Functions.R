#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Testing Thorson Data Import
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.28.17
#
#Purpose: 
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
require(TMB)
require(VAST)
require(FishData)

source('R/load-RACE-data.r')
source('R/create-Data-Geostat.r')

species.codes <- 30060
species.name <- "Sebastes alutus"

working.dir <- getwd()

#============================================================================
# USING FISHSTAT PACKAGE




fd.nr <- download_catch_rates(survey='Aleutian_Islands', add_zeros=TRUE,
                                species_set=species.name)

head(fd.nr)

# fd <- na.omit(fd.nr)
# dim(fd)
# fd.ebs <- download_catch_rates(survey='Eastern_Bering_Sea', add_zeros=FALSE,
#                                species_set='', localdir=paste0(working.dir,'/Data'))
# 
# 
# 
# fd.ai <- download_catch_rates(survey='Aleutian_Islands', add_zeros=FALSE,
#                              species_set=species.codes, localdir=paste0(working.dir,'/Data'))
# 
# fd.goa <- download_catch_rates(survey='Gulf_of_Alaska', add_zeros=FALSE,
#                                species_set=species.codes, localdir=paste0(working.dir,'/Data'))



#============================================================================
# USING ALTERNATIVE FUNCTIONS



c.nr <- create_Data_Geostat(species.codes=species.codes, lat_lon.def='start', survey='AI')

dim(fd.nr)
dim(c.nr)
par(mfrow=c(2,1))
hist(fd.nr$Lat)
hist(c.nr$Lat)

hist(fd.nr$Long)
hist(c.nr$Lon)

#Do some comparative plots
par(mfrow=c(2,1))
hist(fd.nr$Wt, ylab='FishData')
hist(c.nr$Catch_KG)



summary(fd.nr$Wt)
summary(c.nr$Catch_KG)
summary(c.nr$Catch_KG*(0.01/c.nr$AreaSwept_km2))
















