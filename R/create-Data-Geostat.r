
#' Create the Data_Geostat dataframe for vast, based on RACE survey data
#'
#' @param species.codes vector of species codes for which data will be returned
#' @param lat_lon.def string defining how tow-specific Latitude and Longitude will be calculated
#' @param survey string indicating the survey for which data are being extracted: GOA, AI, EBS_SHELF, EBS_SLOPE
#' @param combineSpecies boolean indicating whether species codes should be combined into a single index (i.e. Dusky Rockfish)
#'
#' @return dataframe Data_Geostat with input data for VAST
#' @export
create_Data_Geostat <- function(species.codes, combineSpecies=FALSE, lat_lon.def="start", survey="GOA") {
  ###TESTING###
  # species.codes <- c(30152)#,30420)
  # lat_lon.def <- "mean"
  # survey <- 'EBS_SHELF'
  #############
  
  # source("R/load-RACE-data.r")
  
  #Check Inputs
  if(!lat_lon.def %in% c("mean", "start", "end")) { stop("lat_lon.def must be mean, start, or end") }
  if(!survey %in% c("GOA","AI","EBS_SHELF",'EBS_SLOPE')) { stop(paste("survey is:",survey,", should be one of: GOA, AI, EBS_SHELF, EBS_SLOPE"))  }
  
  #Get Data
  load.data <- load_RACE_data(species.codes=species.codes, combineSpecies=combineSpecies, survey=survey)
  
  #Create VAST input data object
  Data_Geostat <- NULL
  #List Species, if multispecies
  if(length(species.codes) > 1) {
    Data_Geostat$spp <- load.data$Common.Name
  }
  Data_Geostat$Catch_KG <- load.data$WEIGHT  
  Data_Geostat$Year <- load.data$Year 
  Data_Geostat$Vessel <- load.data$VESSEL.x
  Data_Geostat$AreaSwept_km2 <- load.data$effort
  
  #Define Lat and Lon
  if(lat_lon.def=="start") {
    Data_Geostat$Lat <- load.data$START_LATITUDE
    Data_Geostat$Lon <- load.data$START_LONGITUDE
  }
  if(lat_lon.def=="end") {
    Data_Geostat$Lat <- load.data$END_LATITUDE
    Data_Geostat$Lon <- load.data$END_LONGITUDE
  }
  if(lat_lon.def=="mean") {
    Data_Geostat$Lat <- rowMeans(cbind(load.data$START_LATITUDE, load.data$END_LATITUDE), na.rm=TRUE)
    Data_Geostat$Lon <- rowMeans(cbind(load.data$START_LONGITUDE, load.data$END_LONGITUDE), na.rm=TRUE)
  }
  
  #Convert to data frame
  Data_Geostat <- data.frame(Data_Geostat)
  return(Data_Geostat)
}
