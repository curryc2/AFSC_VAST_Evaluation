#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Testing Observation Model Specification
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.30.17
#
#Purpose: To test impact of changing distribution on observation model for positive catch rates,
#           across species.
#  1) Loop through obs model specifications for positive catch rate distribution.
#  2) Run vAST model in parallel across species. 
#
#
#
#==================================================================================================
#NOTES:
#  a) Originally included Normal distribution for positive catch rates,
#       but produced unknown errors for some species.
#==================================================================================================
require(snowfall)
require(parallel)
require(ggplot2)
require(TMB)
require(TMBhelper)
require(VAST)

#Source necessary files
source("R/create-VAST-input.r")
source("R/cleanup-VAST-file.r")

#Create working directory
working.dir <- paste0(getwd(),'/examples/Test_Obs_Model')

#Determine species list
species.list <- read.csv("data/eval_species_list.csv")

#Limit species included
species.list <- species.list[species.list$include=='Y',]
n.species <- nrow(species.list)
#Create
species.series <- c(1:n.species)

#=======================================================================
##### CONTROL SECTION #####
#Number of cores to use
n.cores <- detectCores()-1

#Setup observation model trials
# trial.obs <- c(0:2,5:7)
trial.obs <- c(1:2,5:7)
n.trial.obs <- length(trial.obs)

#Names for trial Obs models
# names.trial.obs <- c('Normal','Lognormal','Gamma','Negative_binomial','Conway-Maxwell_Poisson','Poisson')
names.trial.obs <- c('Lognormal','Gamma','Negative_binomial','Conway-Maxwell_Poisson','Poisson')
#=======================================================================
##### VAST MODEL SPECIFICATIONS #####

Version <- "VAST_v2_4_0"

lat_lon.def <- "mean"

#SPATIAL SETTINGS
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 25
n_x <- c(100, 250, 500, 1000, 2000)[2] # Number of stations
Kmeans_Config <- list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"),
                            west_border = c(-Inf),
                            east_border = c(Inf))


#DERIVED OBJECTS
Region <- "Gulf_of_Alaska"
###########################
# DateFile=paste0(getwd(),'/examples/VAST_output/')

#MODEL SETTINGS
FieldConfig <- c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
RhoConfig <- c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig <- c(Delta1 = 0, Delta2 = 0)


#SPECIFY OUTPUTS
Options <- c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)


#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####


s <- 1
species_wrapper_fxn <- function(s) {
  
  #Define file for analyses
  DateFile <- paste0(trial.dir,"/",species.list$name[s],"/")
  
  #Define species.codes
  species.codes <- species.list$species.code[s]
  
  #=======================================================================
  ##### READ IN DATA AND BUILD VAST INPUT #####
  #  NOTE: this will create the DateFile
  
  VAST_input <- create_VAST_input(species.codes=species.codes, lat_lon.def=lat_lon.def, save.Record=TRUE,
                                  Method="Mesh", grid_size_km=25, n_X=250,
                                  Kmeans_Config=list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 ),
                                  strata.limits=NULL, Region="Gulf_of_Alaska",
                                  DateFile=DateFile,
                                  FieldConfig, RhoConfig, OverdispersionConfig,
                                  ObsModel, Options)
  
  
  
  #Unpack
  TmbData <- VAST_input$TmbData
  Data_Geostat <- VAST_input$Data_Geostat
  Spatial_List <- VAST_input$Spatial_List
  Extrapolation_List <- VAST_input$Extrapolation_List
  
  
  #=======================================================================
  ##### RUN VAST #####
  
  
  
  #Build TMB Object
  #  Compilation may take some time
  TmbList <- VAST::Build_TMB_Fn(TmbData = TmbData, RunDir = DateFile,
                                Version = Version, RhoConfig = RhoConfig, loc_x = Spatial_List$loc_x,
                                Method = Method)
  Obj <- TmbList[["Obj"]]
  
  
  Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                             upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                             bias.correct = FALSE)
  #Save output
  Report = Obj$report()
  Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
  save(Save, file=paste0(DateFile,"Save.RData"))
  
  #========================================================================
  ##### DIAGNOSTIC AND PREDICTION PLOTS #####
  # plot_VAST_output(Opt, Report, DateFile, Region, TmbData, Data_Geostat, Extrapolation_List, Spatial_List)
  
  #========================================================================
  ##### CLEANUP VAST OUTPUT #####
  cleanup_VAST_file(DateFile=DateFile, Version=Version)
  
  #========================================================================
  ##### RETURN SECTION #####
  return(Opt$AIC)
  
} 

#Create a fun output AIC object
AICs <- array(dim=c(n.trial.obs, n.species), dimnames=list(names.trial.obs,species.list$name))


t <- 1
for(t in 1:n.trial.obs) {
  print(paste('## Trial Observation Model',t,'of',n.trial.obs))
  print(paste('# Positive Catch Rate Model is:',names.trial.obs[t]))
  #Specify trial observation model
  ObsModel <- c(trial.obs[t], 0)
  #Setup File
  trial.dir <- paste0(working.dir,"/",names.trial.obs[t])
  dir.create(trial.dir)
  
  #=======================================================================
  ##### SNOWFALL CODE FOR PARALLEL #####
  sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
  sfExportAll() #Exportas all global variables to cores
  sfLibrary(TMB)  #Loads a package on all nodes
  sfLibrary(VAST)
  output <- sfLapply(species.series, fun=species_wrapper_fxn)
  sfStop()
  
  #Output
  AICs[t,] <- unlist(rbind(output))
  
}# next t



