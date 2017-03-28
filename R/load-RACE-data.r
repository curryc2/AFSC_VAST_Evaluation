#' Function to load RACE bottom trawl data. 
#' Catch and haul data are merged with species names, and zeros are added for missing no-catch observations.
#'
#' @param species.codes vector of species codes for which data will be returned
#' @param area string indicating the area of focus: GOA, BS, AI
#' @param writeCSV boolean indicating whether "output/RACE_data_output.csv" is created
#' @param writeDATA boolean indicating whether "output/RACE_data_output.RData" is created
#'
#' @return A data frame of RACE bottom trawl data, with rows equal to species-by-haul observations
#' @export
load_RACE_data <- function(species.codes=c(30152,30420), area="GOA", writeCSV=FALSE, writeDATA=FALSE) {
  require(FishData)

  
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
  catchhaul <- merge(x=catch, y=haul, all.y=TRUE, all.x=FALSE, by.x="HAULJOIN", by.y="HAULJOIN")
  
  #Add in zero observations for catch weight, for no catches.
  #  Drawing on Jim Thorson's code from FishData
  catchhaul.2 <- add_missing_zeros(data_frame=catchhaul, unique_sample_ID_colname="HAULJOIN",
                                   sample_colname="WEIGHT", species_colname="SPECIES_CODE",
                                   species_subset=species.codes,
                                   if_multiple_records="First",
                                   Method="Fast")
  #Limit to specific area
  catchhaul.2 <- catchhaul.2[catchhaul.2$REGION.x==area,]
  
  #Calculate and add Effort and CPUE
  catchhaul.2$effort <- catchhaul.2$NET_WIDTH*catchhaul.2$DISTANCE_FISHED/1000
  catchhaul.2$cpue <- catchhaul.2$WEIGHT/catchhaul.2$effort
  
  #Add species name
  species.code.data <- read.csv("data/race_species_codes.csv", header=TRUE, stringsAsFactors=FALSE)
  catchhaul.3 <- merge(x=catchhaul.2, y=species.code.data[,c("Species.Code","Common.Name")], 
                       by.x="SPECIES_CODE", by.y="Species.Code")
  
  #Add year of survey
  cruise.info <- read.csv("data/race_cruise_info.csv", header=TRUE, stringsAsFactors=FALSE)
  output <- merge(x=catchhaul.3, y=cruise.info[,c("Cruise.Join.ID","Year")], 
                       by.x="CRUISEJOIN.x", by.y="Cruise.Join.ID")
  
  #Return Section
  if(writeCSV==TRUE) { write.csv(output, file="output/RACE_data_output.csv") }
  if(writeDATA==TRUE) { save(output, file="output/RACE_data_output.RData") }
  
  return(output)
}