#######
#' Create input objects for VAST model
#'   NOTE: Currently only tested for single-species application
#'
#' @param species.codes vector of species codes to be evaluated
#' @param lat_lon.def string defining how tow-specific Latitude and Longitude will be calculated
#' @param Method 
#' @param grid_size_km 
#' @param n_x 
#' @param Kmeans_Config 
#' @param strata.limits dataframe of strata limits for post-hoc apportionment
#' @param survey string indicating the survey for which data are being extracted: GOA, AI, EBS_SHELF, EBS_SLOPE
#' @param DateFile path for directory housing VAST model and output figures and objects
#' @param FieldConfig 
#' @param RhoConfig 
#' @param OverdispersionConfig 
#' @param ObsModel 
#' @param Options 
#' @param save.Record boolean indicating whether or not VAST settings record is saved
#'
#' @return VAST_input: Containing Data_Geostat, Spatial_List, Extrapolation_List, and TmbData
#' @export
create_VAST_input <- function(species.codes, combineSpecies=FALSE, lat_lon.def="mean", save.Record=TRUE,
                                Method="Mesh", grid_size_km=25, n_x=250,
                                Kmeans_Config=list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 ),
                                strata.limits=NULL, survey="GOA",
                                DateFile=paste0(getwd(),'/VAST_output/'),
                                FieldConfig, RhoConfig, OverdispersionConfig,
                                ObsModel, Options) {
  
  # source('R/create-Data-Geostat.r')
  
  #vERSION NUMBER
  Version  <- "VAST_v2_4_0"
  
  #DATA Set
  Data_Set <- "VAST_EVAL"
  
  #Operation if no strata.limits are defined
  if(is.null(strata.limits)) {
    strata.limits <- data.frame(STRATA = c("All_areas"))#,
                                # west_border = c(-Inf),
                                # east_border = c(Inf))
  }
  
  #Determine Correct area to allign with RACE survey
  if(survey %in% c("GOA","AI","EBS_SHELF",'EBS_SLOPE')) { 
    if(survey=="GOA") { Region <- "Gulf_of_Alaska"; area <- "GOA" }
    if(survey=="AI") { Region <- "Aleutian_Islands"; area <- "AI" }
    if(survey=="EBS_SHELF" | survey=="EBS_SLOPE") { Region <- "Eastern_Bering_Sea"; area <- "BS" }
    
  }else {
    stop(paste("survey is:",survey,", should be one of: GOA, AI, EBS_SHELF, EBS_SLOPE"))
  }
  
  # Region <- "Other"
  
  #Retreive Data
  Data_Geostat <- create_Data_Geostat(species.codes=species.codes, lat_lon.def=lat_lon.def, survey=survey) 
  
  #Build Extrapolition Grid
  start.time <- date()
  if(Region=="Other") {
    Extrapolation_List  <- SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region=Region, strata.limits=strata.limits,
                                                                             observations_LL=Data_Geostat[,c("Lat","Lon")],
                                                                             maximum_distance_from_sample=15)
  }else {
    Extrapolation_List  <- SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region=Region, strata.limits=strata.limits)
  }
  end.time <- date()
  #Create Location for Saving Files
  dir.create(DateFile) #Recursive may need to be false if other elements exist
  
  #Save Settings for Later reference
  if(save.Record==TRUE) {
    warning("Currently save.Record=FALSE only accepted, easy-read records will not be saved.")
    # Record = ThorsonUtilities::bundlelist(c("Data_Set"=Data_Set,
    #                                         "Version"=Version, "Method"=Method, "grid_size_km"=grid_size_km,
    #                                         "n_x"=n_X, "FieldConfig"=FieldConfig,
    #                                         "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig,
    #                                         "ObsModel"=ObsModel,
    #                                         "Kmeans_Config"=Kmeans_Config))
    # # save(Record, file = file.path(DateFile, "Record.RData"))
    # capture.output(Record, file = paste0(DateFile, "Record.txt"))
  
  }
  #Generate Information used for spatio-temporal estimation
  Spatial_List <- SpatialDeltaGLMM::Spatial_Information_Fn(grid_size_km = grid_size_km,
                                         n_x = n_x, Method = Method, Lon = Data_Geostat[,"Lon"], Lat = Data_Geostat[, "Lat"],
                                         Extrapolation_List = Extrapolation_List,
                                         randomseed = Kmeans_Config[["randomseed"]],
                                         nstart = Kmeans_Config[["nstart"]],
                                         iter.max = Kmeans_Config[["iter.max"]], DirPath = DateFile,
                                         Save_Results = FALSE)
  #Update Data_Geostat
  Data_Geostat <- cbind(Data_Geostat, knot_i = Spatial_List$knot_i)
  
  #Build VAST model data input
  
  if(length(species.codes) > 1) {
    #MULTISPECIES
    TmbData <- VAST::Data_Fn(Version=Version, FieldConfig=FieldConfig, 
                      OverdispersionConfig=OverdispersionConfig,
                      RhoConfig=RhoConfig, ObsModel=ObsModel, c_i=as.numeric(Data_Geostat[,'spp'])-1,
                      b_i=Data_Geostat[,'Catch_KG'], a_i=Data_Geostat[,'AreaSwept_km2'],
                      v_i=as.numeric(Data_Geostat[,'Vessel'])-1, s_i=Data_Geostat[,'knot_i']-1,
                      t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl, MeshList=Spatial_List$MeshList,
                      GridList=Spatial_List$GridList, Method=Spatial_List$Method, Options=Options )
  }else {
    #SINGLE SPECIES
    TmbData = VAST::Data_Fn(Version = Version, FieldConfig = FieldConfig,
                      OverdispersionConfig = OverdispersionConfig, RhoConfig = RhoConfig,
                      ObsModel = ObsModel, c_i = rep(0, nrow(Data_Geostat)),
                      b_i = Data_Geostat[, "Catch_KG"], a_i = Data_Geostat[,"AreaSwept_km2"],
                      v_i = as.numeric(Data_Geostat[,"Vessel"]) - 1,
                      s_i = Data_Geostat[, "knot_i"] - 1, t_i = Data_Geostat[, "Year"],
                      a_xl = Spatial_List$a_xl,
                      MeshList = Spatial_List$MeshList, GridList = Spatial_List$GridList,
                      Method = Spatial_List$Method, Options = Options)
  }

  
  
  #Return Section
  VAST_input <- NULL
  VAST_input$Data_Geostat <- Data_Geostat
  VAST_input$Spatial_List <- Spatial_List
  VAST_input$TmbData <- TmbData
  VAST_input$Extrapolation_List <- Extrapolation_List
  
  return(VAST_input)
}



