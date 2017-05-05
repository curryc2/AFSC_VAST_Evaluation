#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Knot Number
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 4.10.17
#
#Purpose: To explore sensitivity of model-based index estimates to different knot number specification
#
#
#
#==================================================================================================
#NOTES:
#
#  a) Could calculate design-based estiamtes within wrapper function or at the end while plotting
#  
#==================================================================================================
#TIMING:
# For 10-1000, by 100 knots
# [1] "### START: Mon Apr 10 16:17:08 2017"
# [1] "### END: Mon Apr 10 22:54:04 2017"
#==================================================================================================

require(VAST)
require(TMB)
require(parallel)
require(snowfall)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/cleanup-VAST-file.r")



#Create working directory
working.dir <- paste0(getwd(),"/examples/Test_Knot_Number")

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

#Boolean for running estimation models
do.estim <- FALSE

#Trial Knot Numbers
trial.knots <- seq(100, 1000, by=100)
n.trial.knots <- length(trial.knots)

#Boolean for bias correction
bias.correct <- FALSE
#=======================================================================
##### Run VAST model  #####
Version <- "VAST_v2_4_0"
Region <- "Gulf_of_Alaska"
lat_lon.def <- "mean"

#SPATIAL SETTINGS
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
# n_x = c(100, 250, 500, 1000, 2000)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"),
                            west_border = c(-Inf),
                            east_border = c(Inf))

#MODEL SETTINGS
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)



#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####


s <- 1
# for(s in 1:n.species) {
species_wrapper_fxn_knots <- function(s, n_x) {
  
  #Define file for analyses
  DateFile <- paste0(trial.dir,"/",species.list$name[s],"/")
  dir.create(DateFile)
  #Define species.codes
  species.codes <- species.list$species.code[s]
  
  #Calculate design-based estimate
  db_est <- calc_design_based_index(species.codes=species.codes, Region=Region)
  save(db_est, file=paste0(DateFile,"db_est.RData"))
  #=======================================================================
  ##### READ IN DATA AND BUILD VAST INPUT #####
  #  NOTE: this will create the DateFile
  
  VAST_input <- create_VAST_input(species.codes=species.codes, lat_lon.def=lat_lon.def, save.Record=FALSE,
                                  Method=Method, grid_size_km=grid_size_km, n_x=n_x,
                                  Kmeans_Config=Kmeans_Config,
                                  strata.limits=strata.limits, Region=Region,
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
                             bias.correct = bias.correct)
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
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List",
     "TmbList", "Obj", "Opt", "Report", "Save",
     "db_est")
  
  #========================================================================
  ##### RETURN SECTION #####
  # return(Opt$AIC)
  
} 


#=======================================================================
##### Loop Through Trial Knots  #####
if(do.estim==TRUE) {
  time.1 <- date()
  
  t <- 1
  for(t in 1:n.trial.knots) {
    print(paste('## Trial Knot Number',t,'of',n.trial.knots))
    print(paste('# Trial Knots:',trial.knots[t]))
    #Specify trial observation model

    #Specify knots
    n_x <- trial.knots[t]
    
    #Setup File
    trial.dir <- paste0(working.dir,"/",n_x,"_bias.corr_",bias.correct)
    dir.create(trial.dir)
    

    
    #=======================================================================
    ##### SNOWFALL CODE FOR PARALLEL #####
    sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
    sfExportAll() #Exportas all global variables to cores
    sfLibrary(TMB)  #Loads a package on all nodes
    sfLibrary(VAST)
    output <- sfLapply(species.series, fun=species_wrapper_fxn_knots, n_x=n_x)
    # sfRemove(Save)
    # sfRemover(VAST_input)
    sfStop()
    
    
    
  }# next t



  #Go back through and delete Unnecessary parameter estimates
  t <- 1
  for(t in 1:n.trial.knots) {
    s <- 1
    for(s in 1:n.species) {
      temp.DateFile <- paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct,"/",species.list$name[s],"/")
      #Remove manually
      file.remove(paste0(temp.DateFile,"parameter_estimates.Rdata"))
      file.remove(paste0(temp.DateFile,"parameter_estimates.txt"))
    }#next s
  }#next t
  
  time.2 <- date()

  print(paste('### START:', time.1))
  print(paste('### END:', time.2))

}
#=======================================================================
##### Plot Comparison of Results #####=
#Flag for excluding the 2001 design-based index from comparison,
#  given unequal spatial sampling.
exclude.2001 <- TRUE 

n.species
n.trial.knots


#Create data objects
#First get years
load(paste0(working.dir,"/",trial.knots[1],"_bias.corr_",bias.correct,"/",species.list$name[1],"/db_est.RData"))
yrs.surv <- db_est$YEAR
n.yrs.surv <- length(yrs.surv)

#Read in data


#Plot Results









