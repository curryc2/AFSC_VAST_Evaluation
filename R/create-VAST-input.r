create_VAST_input <- function(species.codes, lat_lon.def="mean") {
  #TESTING
  species.codes <- c(30152,30420)
  lat_lon.def <- "mean"
  
  source("R/load-RACE-data.r")
  
  #Check Inputs
  if(!lat_lon.def %in% c("mean", "start", "end")) { stop("lat_lon.def must be mean, start, or end") }
  
  #Get Data
  load.data <- load_RACE_data(species.codes=species.codes, area="GOA")
  
  #Create VAST input data object
  Data_Geostat <- NULL
  Data_Geostat$Catch_KG <- temp.data$WEIGHT  
  Data_Geostat$Year <- year(temp.data$START_TIME) 
}