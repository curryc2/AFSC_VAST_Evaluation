#' Function to load RACE bottom trawl data. 
#' Catch and haul data are merged with species names, and zeros are added for missing no-catch observations.
#'
#' @param species.codes vector of species codes for which data will be returned
#' @param survey string indicating the survey for which data are being extracted: GOA, AI, EBS_SHELF, EBS_SLOPE
#' @param writeCSV boolean indicating whether "output/RACE_data_output.csv" is created
#' @param writeDATA boolean indicating whether "output/RACE_data_output.RData" is created
#'
#' @return A data frame of RACE bottom trawl data, with rows equal to species-by-haul observations
#' @export
load_RACE_data <- function(species.codes=c(30152,30420), survey="GOA", writeCSV=FALSE, writeDATA=FALSE) {
  require(FishData)
  require(dplyr)
  ###TESTING###
  species.codes <- 30420#21720 #Pacific Cod 
  survey <- 'GOA'
  #############
  
  if(!survey %in% c("GOA","AI","EBS_SHELF",'EBS_SLOPE')) { stop(paste("survey is:",survey,"should be one of: GOA, AI, EBS_SHELF, EBS_SLOPE"))  }
  #FOR ORIGINAL .CSV INPUT
  # #Opening section to determine suitability
  # if(file.exists("data/race_base_haul.csv")==FALSE) { stop("data/race_base_haul.csv NOT FOUND") }
  # if(file.exists("data/race_base_catch.csv")==FALSE) { stop("data/race_base_catch.csv NOT FOUND") }
  # if(file.exists("data/race_species_codes.csv")==FALSE) { stop("data/race_species_codes.csv NOT FOUND") }
  # 
  # #Load catch and haul data
  # haul <- read.csv("data/race_base_haul.csv", header=TRUE, stringsAsFactors=FALSE)
  # catch <- read.csv("data/race_base_catch.csv", header=TRUE, stringsAsFactors=FALSE)
  
  #FOR NEW .RData INPUT
  #Opening section to determine suitability
  if(file.exists("data/race_base_haul.RData")==FALSE) { stop("data/race_base_haul.RData NOT FOUND") }
  if(file.exists("data/race_base_catch.RData")==FALSE) { stop("data/race_base_catch.RData NOT FOUND") }
  if(file.exists("data/race_species_codes.csv")==FALSE) { stop("data/race_species_codes.csv NOT FOUND") }
  if(file.exists("data/race_cruise_info.csv")==FALSE) { stop("data/race_cruise_info.csv NOT FOUND") }
  
  #Load catch and haul data
  load("data/race_base_catch.RData")
  load("data/race_base_haul.RData")

  
  #Merge data
  # catchhaul <- merge(x=catch, y=haul, all.y=TRUE, all.x=FALSE, by.x="HAULJOIN", by.y="HAULJOIN")
  #dplyr it for greater speed
  catchhaul <- right_join(x=catch, y=haul, by=c("HAULJOIN"))
  
  #Add in zero observations for catch weight, for no catches.
  #  Drawing on Jim Thorson's code from FishData
  catchhaul.2 <- FishData::add_missing_zeros(data_frame=catchhaul, unique_sample_ID_colname="HAULJOIN",
                                               sample_colname="WEIGHT", species_colname="SPECIES_CODE",
                                               species_subset=species.codes,
                                               if_multiple_records="First",
                                                Method="Fast")
  dim(catchhaul.2)
  
  #Bring in cruise info
  #Add year of survey and name
  cruise.info <- read.csv("data/race_cruise_info.csv", header=TRUE, stringsAsFactors=FALSE)
  catchhaul.3 <- left_join(x=catchhaul.2, y=cruise.info[,c("Cruise.Join.ID","Year","Survey")],
                             by=c("CRUISEJOIN.x"="Cruise.Join.ID"))
  
  
  
  # #Limit to specific survey
  # catchhaul.3 <- catchhaul.3[catchhaul.3$Survey==survey,]
  
  #Calculate and add Effort and CPUE
  catchhaul.3$effort <- catchhaul.3$NET_WIDTH*catchhaul.3$DISTANCE_FISHED/1000
  catchhaul.3$cpue <- catchhaul.3$WEIGHT/catchhaul.3$effort
  
  #Add species name
  species.code.data <- read.csv("data/race_species_codes.csv", header=TRUE, stringsAsFactors=FALSE)
  output <- left_join(x=catchhaul.3, y=species.code.data[,c("Species.Code","Common.Name")], 
                       by=c("SPECIES_CODE"="Species.Code"))
  output.2 <- merge(x=catchhaul.3, y=species.code.data[,c("Species.Code","Common.Name")], 
                       by.x="SPECIES_CODE", by.y="Species.Code")
  
  # #Limit to specific survey
  catchhaul.3 <- catchhaul.3[catchhaul.3$Survey==survey,]
  
  #Return Section
  if(writeCSV==TRUE) { write.csv(output, file="output/RACE_data_output.csv") }
  if(writeDATA==TRUE) { save(output, file="output/RACE_data_output.RData") }
  
  return(output)
}

