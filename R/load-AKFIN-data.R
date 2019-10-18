#' Function to load RACE bottom trawl data. 
#' Catch and haul data are merged with species names, and zeros are added for missing no-catch observations.
#'
#' @param species.codes vector of species codes for which data will be returned
#' @param survey string indicating the survey for which data are being extracted: GOA, AI, EBS_SHELF, EBS_SLOPE
#' @param writeCSV boolean indicating whether "output/RACE_data_output.csv" is created
#' @param writeDATA boolean indicating whether "output/RACE_data_output.RData" is created
#' @param combineSpecies boolean indicating whether species codes should be combined into a single index (i.e. Dusky Rockfish)
#'
#' @return A data frame of RACE bottom trawl data, with rows equal to species-by-haul observations
#' @export
load_AKFIN_data <- function(species.codes=c(30150,30152), combineSpecies=FALSE, survey="GOA", writeCSV=FALSE, writeDATA=FALSE) {
  require(FishData)
  require(dplyr)
  require(readr)
  ##TESTING###
  species.codes <- 30420
  combineSpecies <- FALSE
  survey <- 'GOA'
  writeCSV <- FALSE
  writeDATA <- FALSE
  #############
  
  if(survey %in% c("GOA","AI","EBS_SHELF",'EBS_SLOPE')) { 
    if(survey=="GOA") { Region <- "Gulf_of_Alaska"; area <- "GOA" }
    if(survey=="AI") { Region <- "Aleutian_Islands"; area <- "AI" }
    if(survey=="EBS_SHELF" | survey=="EBS_SLOPE") { Region <- "Eastern_Bering_Sea"; area <- "BS" }
    
  }else {
    stop(paste("survey is:",survey,", should be one of: GOA, AI, EBS_SHELF, EBS_SLOPE"))
  }
  #FOR ORIGINAL .CSV INPUT
  akfin.haul <- readr::read_csv(file.path("data","AKFIN Fall 2019","Haul Descriptions.csv"))
  head(akfin.haul)
  akfin.cpue <- readr::read_csv(file.path("data","AKFIN Fall 2019","race_cpue_by_haul.csv"), skip=5)
  head(akfin.cpue)
  
  # Save as .rds for quicker read-in
  # saveRDS(akfin.haul, )
  
  # akfin.catch <- readr::read_csv(file.path("data","AKFIN Fall 2019","race_catch_by_haul.csv"), skip=5)
  # head(akfin.catch)
  
  
    #FOR NEW .RData INPUT
    #Opening section to determine suitability
    # if(file.exists("data/race_base_haul.RData")==FALSE) { stop("data/race_base_haul.RData NOT FOUND") }
    # if(file.exists("data/race_base_catch.RData")==FALSE) { stop("data/race_base_catch.RData NOT FOUND") }
    if(file.exists("data/race_base_haul.rds")==FALSE) { stop("data/race_base_haul.rds NOT FOUND") }
    if(file.exists("data/race_base_catch.rds")==FALSE) { stop("data/race_base_catch.rds NOT FOUND") }
    if(file.exists("data/race_species_codes.csv")==FALSE) { stop("data/race_species_codes.csv NOT FOUND") }
    if(file.exists("data/race_cruise_info.csv")==FALSE) { stop("data/race_cruise_info.csv NOT FOUND") }
    
    #Load catch and haul data
    # load("data/race_base_catch.RData")
    # load("data/race_base_haul.RData")
    catch <- readRDS("data/race_base_catch.rds")
    haul <- readRDS("data/race_base_haul.rds")
    #Limit to only certified abundance hauls
    haul <- haul[haul$ABUNDANCE_HAUL=='Y',]
    
    
    #Merge data
    # catchhaul <- merge(x=catch, y=haul, all.y=TRUE, all.x=FALSE, by.x="HAULJOIN", by.y="HAULJOIN")
    #dplyr it for greater speed
    catchhaul <- right_join(x=catch, y=haul, by=c("HAULJOIN"))
    
    #Add in zero observations for catch weight, for no catches.
    #  Drawing on Jim Thorson's code from FishData
    #NOTE: species.codes are now treated as a factor.
    catchhaul.2 <- FishData::add_missing_zeros(data_frame=catchhaul, unique_sample_ID_colname="HAULJOIN",
                                               sample_colname="WEIGHT", species_colname="SPECIES_CODE",
                                               species_subset=species.codes,
                                               if_multiple_records="First",
                                               Method="Fast")
    
    
    
    
    #Bring in cruise info
    #Add year of survey and name
    cruise.info <- read.csv("data/race_cruise_info.csv", header=TRUE, stringsAsFactors=FALSE)
    # catchhaul.3 <- left_join(x=catchhaul.2, y=cruise.info[,c("Cruise.Join.ID","Year","Survey")],
    #                            by=c("CRUISEJOIN.x"="Cruise.Join.ID"))
    
    #User inner_join because it removes hauls WITHOUT identified Year and Survey
    catchhaul.3 <- inner_join(x=catchhaul.2, y=cruise.info[,c("Cruise.Join.ID","Year","Survey")],
                              by=c("CRUISEJOIN.x"="Cruise.Join.ID"))
    
    
    # #Find difference
    # loc <- which(catchhaul.3$CRUISEJOIN.x %in% catchhaul.3.in$CRUISEJOIN.x)
    # catchhaul.3[-loc,]
  }else {
    
    
    # Using AKFIN DATA ===========================================
    akfin.haul <- readr::read_csv(file.path("data","AKFIN Fall 2019","Haul Descriptions.csv"))
    head(akfin.haul)
    akfin.cpue <- readr::read_csv(file.path("data","AKFIN Fall 2019","race_cpue_by_haul.csv"), skip=5)
    head(akfin.cpue)
    # akfin.catch <- readr::read_csv(file.path("data","AKFIN Fall 2019","race_catch_by_haul.csv"), skip=5)
    # head(akfin.catch)
    
    akfin.haul <- akfin.haul
    
    #=============================================================
    
  }
  #Limit to specific survey
  catchhaul.3 <- catchhaul.3[catchhaul.3$Survey==survey,]
  # catchhaul.3 <- catchhaul.3[catchhaul.3$Survey==survey & catchhaul.3$REGION.x==area,]
  
  #=================================================
  #AGGREGATE CATCH BIOMASS ACROSS SPECIES CODES IN THE CASE OF A COMBINED INDEX
  #  Dusky Rockfish Example 2 species codes to single index
  #   Variables to combine: WEIGHT
  if(combineSpecies==TRUE) {
    catchhaul.4 <- data.frame(catchhaul.3 %>% group_by(HAULJOIN) %>% 
                                mutate('WEIGHT'=sum(WEIGHT, na.rm=TRUE)))
    #Since we have aggregated, only retain rows for 1st listed species code
    catchhaul.5 <- catchhaul.4[catchhaul.4$SPECIES_CODE==species.codes[1],]
  }else {
    catchhaul.5 <- catchhaul.3
  }
  #=================================================
  
  
  #Calculate and add Effort and CPUE
  catchhaul.5$effort <- catchhaul.5$NET_WIDTH*catchhaul.5$DISTANCE_FISHED/1000
  catchhaul.5$cpue <- catchhaul.5$WEIGHT/catchhaul.5$effort
  
  #Add species name
  species.code.data <- read.csv("data/race_species_codes.csv", header=TRUE, stringsAsFactors=FALSE)
  
  #Convert codes to a factor
  # species.code.data$Species.Code <- as.factor(species.code.data$Species.Code)
  # 
  # output <- left_join(x=catchhaul.3, y=species.code.data[,c("Species.Code","Common.Name")],
  #                      by=c("SPECIES_CODE"="Species.Code"))
  output <- merge(x=catchhaul.5, y=species.code.data[,c("Species.Code","Common.Name")], 
                  by.x="SPECIES_CODE", by.y="Species.Code")
  
  #Return Section
  if(writeCSV==TRUE) { write.csv(output, file="output/RACE_data_output.csv") }
  if(writeDATA==TRUE) { save(output, file="output/RACE_data_output.RData") }
  
  return(output)
}


##### TESTING FUNCTION #####
# #Northern Rockfish
# temp <- load_RACE_data(species.codes=c(30420), combineSpecies=FALSE, survey='GOA', writeCSV=FALSE, writeDATA=FALSE)
# dim(temp)
# 
# #Dusky Rockfish partial
# temp.2 <- load_RACE_data(species.codes=c(30150), combineSpecies=FALSE, survey='GOA', writeCSV=FALSE, writeDATA=FALSE) 
# dim(temp.2)
# 
# 
# #Dusky Rockfish Complete/combined
# temp.3 <- load_RACE_data(species.codes=c(30150,30152), combineSpecies=TRUE, survey='GOA', writeCSV=FALSE, writeDATA=FALSE)
# dim(temp.3)
# 
# #Checkup
# unique(sort(unique(temp.3$HAULJOIN))==sort(unique(temp$HAULJOIN)))
# unique(sort(unique(temp.2$HAULJOIN))==sort(unique(temp$HAULJOIN)))
# 
# 
# #Checkup figure
# require(ggplot2)
# require(ggthemes)
# 
# temp.4 <- load_RACE_data(species.codes=c(30150,30152), combineSpecies=FALSE, survey='GOA', writeCSV=FALSE, writeDATA=FALSE)
# dim(temp.4)
# 
# 
# 
# g <- ggplot(temp.3, aes(x=Year, y=WEIGHT)) +
#        stat_summary(fun.y=sum, geom='area', fill='blue') +
#        stat_summary(fun.y=sum, geom='point', color='red')
# 
# g
# 
# g2 <- ggplot(temp.4, aes(x=Year, y=WEIGHT, fill=SPECIES_CODE)) +
#         stat_summary(fun.y=sum, geom='area') +
#         stat_summary(fun.y=sum, geom='point', color='black') +
#         facet_wrap(~SPECIES_CODE, ncol=1)
# 
# g2
# 
# g3 <- ggplot(temp.4, aes(x=Year, y=WEIGHT, fill=SPECIES_CODE)) +
#         stat_summary(fun.y=sum, geom='area', aes(group=SPECIES_CODE), alpha=0.5) +
#         stat_summary(fun.y=sum, geom='point', color='black') +
#         scale_fill_colorblind()
# 
# g3




