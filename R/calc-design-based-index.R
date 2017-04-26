#' Function to calculate design-based survey abundance estimate from RACE survey data.
#'
#' @param species.codes vector of species codes for which data will be returned
#' @param Region string indicating region of evaluation, one of: Gulf_of_Alaska, Eastern_Bering_Sea, or Aleutian_Islands
#'
#' @return data frame containing annual survey biomass estimate, variance, SD, and CV.
#' @export
calc_design_based_index <- function(species.codes, Region="Gulf_of_Alaska") {
  ### TESTING VALUES
  species.codes <- 21720 #Pacific Cod
  Region <- 'Gulf_of_Alaska'
  ###
  require(dplyr)
  source("R/load-RACE-data.r")
  #Input Checking...
  if(Region %in% c("Gulf_of_Alaska", "Eastern_Bering_Sea", "Aleutian_Islands")) {
    if(Region=="Gulf_of_Alaska") { area <- "GOA"; survey <- "GOA" }
    if(Region=="Eastern_Bering_Sea") { area <- "BS"; stop("Error: Currently not implemented for EBS, need to correct slope vs shelf survey") }
    if(Region=="Aleutian_Islands") { area <- "AI"; survey <- "AI" }
  }else {
    stop("Region must be one of: Gulf_of_Alaska, Eastern_Bering_Sea, Aleutian_Islands")
  }
  if(length(species.codes)>1) { stop("Error: calc_design_based_index is currently only tested for a single species.") }
  
  #Load RACE survey data
  load.data <- load_RACE_data(species.codes=species.codes, area=area, writeCSV=FALSE, writeDATA=FALSE)
  
  #Calculate design-based estimator
  
  strata <- sort(unique(load.data$STRATUM))
  n.strata <- length(strata)
  
  #Calculate sum and var in cpue by year and stratum
  # cstrat<-ddply(load.data ,c("Year","STRATUM"), summarize, CPUE=sum(cpue), CPUEvar=var(cpue)) #plyr
  cstrat <- data.frame(load.data %>% group_by(Year, STRATUM) %>% summarize(CPUE=sum(cpue), CPUEvar=var(cpue))) #dplyr
  
  #Calculate number of hauls in each year and stratum
  # hstrat <- ddply(load.data, c("Year","STRATUM"), summarize, n_sta=length(unique(HAULJOIN))) #plyr
  hstrat <- data.frame(load.data %>% group_by(Year, STRATUM) %>% summarize(n_sta=length(unique(HAULJOIN)))) #dplyr
  
  #Join together
  # biomvar <- merge(cstrat, hstrat, by.x=c("Year","STRATUM"), by.y=c("Year","STRATUM"), all.x=TRUE) #all.x=TRUE, all.y=FALSE
  biomvar <- left_join(cstrat, hstrat, by=c("Year", "STRATUM")) #dplyr
  colnames(biomvar)<-c("YEAR","STRATUM","CPUE","VAR","n_stations")
  
  
  #Load Strata Data
  strata.data <- read.csv("data/race_stratum_info.csv", header=TRUE)
  #Limit to correct survey area
  strata.area <- strata.data[strata.data$Survey==area,c(2,3,5)] #NEEDS TO BE UPDATED
  names(strata.area) <- c("STRATUM","AREA","INPFC_AREA")
  
  # biomvar <- merge(biomvar, strata.area, by=c("STRATUM"), all.x=TRUE) #plyr
  biomvar <- left_join(biomvar, strata.area, by=c("STRATUM")) #dplyr
  
  biomvar$BIOMASS<-(biomvar$CPUE/biomvar$n_stations)*biomvar$AREA
  biomvar$VAR2<-biomvar$AREA^2*(biomvar$VAR/biomvar$n_stations)
  
  #Calculate complete
  biomass <- data.frame(biomvar %>% group_by(YEAR) %>%
                          summarize(Biomass=sum(BIOMASS,na.rm=TRUE)/1e3,
                                    Variance=sum(VAR2,na.rm=TRUE)/(1e3^2),
                                    SD=sqrt(sum(VAR2,na.rm=TRUE))/1e3,
                                    CV=SD/(sum(BIOMASS,na.rm=TRUE)/1e3)))
  
  return(biomass)
  
}