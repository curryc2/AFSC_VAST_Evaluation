#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Temporal Linkages Between
#                                                              Intercepts and Spatio-temporal RE
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 4.10.17
#
#Purpose: To explore VAST model-based index sensitivity to autoregressive specifications on:
#           1) Intercept
#           2) Spatio-temporal random effect
#           3) Fix the number of knots @ 500, to begin with and consider changing to c(100,500,1000) for sensitivity
#
#
#
#==================================================================================================
#NOTES:
#
#  
#==================================================================================================
#TIMING:
#
#==================================================================================================
require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Autoregressive")


#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

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
do.estim <- TRUE

#Trial Knot Numbers
trial.knots <- c(100,500,1000)
n.trial.knots <- length(trial.knots)

#Trial AUTOREGRESSIVE specifications
#Note starts at 0
# rho.int.types <- c('Fixed_Effect','Independent_Among_Years','Random_Walk','Constant_Intercept','Autoregressive')
rho.int.types <- c('FE','IaY','RW','CI','AR')
# rho.stRE.types <- c('Independent_Among_Years',NA,'Random_Walk',NA,'Autoregressive')
rho.stRE.types <- c('IaY',NA,'RW',NA,'AR')

#Read in Autoregressive Input
trial.rho <- t(read.csv('Data/Test-Autoregressive-Input.csv', header=TRUE, stringsAsFactors=FALSE)[,-c(1:2)])
n.trial.rho <- nrow(trial.rho)

# #Intercept
# rho.inter.ep <- c(0) #Encounter Probability Component 
# rho.inter.pcr <- c(0) #Positive Catch Rate Component
# #Spatio-temporal RE
# rho.stRE.ep #Encounter Probability Component 
# rho.stRE.pcr #Positive Catch Rate Component

#Boolean for bias correction
bias.correct <- FALSE
#=======================================================================
##### Run VAST model  #####
Version <- "VAST_v2_4_0"
lat_lon.def <- "start"

#SPATIAL SETTINGS
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
# n_x = c(100, 250, 500, 1000, 2000)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"))

#MODEL SETTINGS
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
# RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#Output Directory Name
output.dir <- paste0(working.dir,"/output_bias.correct_",bias.correct)



#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####

s <- 1
# for(s in 1:n.species) {
wrapper_fxn <- function(s, n_x, RhoConfig) {
  
  #Define file for analyses
  DateFile <- paste0(trial.dir,"/",species.list$survey[s],"_",species.list$name[s],"/")
  
  dir.create(DateFile)
  
  #Define species.codes
  species.codes <- species.list$species.code[s]
  survey <- species.list$survey[s]
  
  #=======================================================================
  ##### READ IN DATA AND BUILD VAST INPUT #####
  #  NOTE: this will create the DateFile
  
  VAST_input <- create_VAST_input(species.codes=species.codes, lat_lon.def=lat_lon.def, save.Record=FALSE,
                                  Method=Method, grid_size_km=grid_size_km, n_x=n_x,
                                  Kmeans_Config=Kmeans_Config,
                                  strata.limits=strata.limits, survey=survey,
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
  # Report = Obj$report()
  Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
  # save(Save, file=paste0(DateFile,"Save.RData"))
  
  #Calculate index values
  # TmbData = TmbData, Sdreport = Opt[["SD"]]
  vast_est <- get_VAST_index(TmbData=TmbData, Sdreport=Opt[["SD"]], bias.correct=bias.correct, Data_Geostat=Data_Geostat)
  #========================================================================
  ##### DIAGNOSTIC AND PREDICTION PLOTS #####
  # plot_VAST_output(Opt, Report, DateFile, survey, TmbData, Data_Geostat, Extrapolation_List, Spatial_List)
  
  #========================================================================
  ##### CLEANUP VAST OUTPUT #####
  # cleanup_VAST_file(DateFile=DateFile, Version=Version) #No longer necessary as we are deleting everything at the end
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List",
     "TmbList", "Obj", "Save")#, "Opt", "Report")
  
  #========================================================================
  setwd(home.dir)
  ##### RETURN SECTION #####
  out <- NULL
  out$vast_est <- vast_est
  out$Opt <- Opt
  out$Report <- Report
  return(out)
} 


#=======================================================================
##### Loop Through Trial Knots  #####
if(do.estim==TRUE) {
  vast_est.output <- vector('list', length=(n.trial.knots * n.trial.rho))
  vast_knots <- vector(length=(n.trial.knots * n.trial.rho))
  vast_rho.int <- vector(length=(n.trial.knots * n.trial.rho))
  vast_rho.stRE <- vector(length=(n.trial.knots * n.trial.rho))
  
  time.1 <- date()
  
  #Counter for knots by rho
  counter <- 1
  
  t <- 1
  for(t in 1:n.trial.knots) {
    print(paste('## Trial Knot Number',t,'of',n.trial.knots))
    print(paste('# Trial Knots:',trial.knots[t]))
    #Specify trial observation model
    
    #Specify knots
    n_x <- trial.knots[t]
    
    r <- 1
    for(r in 1:n.trial.rho) {
      #Specify intercepts and spatio-temporal variation across time
      RhoConfig <- trial.rho[r,]
      names(RhoConfig) <- c('Beta1','Beta2','Epsilon1','Epsilon2')
      
      #Record
      if(RhoConfig[1]==RhoConfig[2]) {#IF intercept specs are the same
        vast_rho.int[counter] <- rho.int.types[RhoConfig[1]+1]
      }else {#IF different
        vast_rho.int[counter] <- paste0(rho.int.types[RhoConfig[1]+1], "-", rho.int.types[RhoConfig[2]+1])
      }
      if(RhoConfig[3]==RhoConfig[4]) {
        vast_rho.stRE[counter] <- rho.stRE.types[RhoConfig[3]+1]
      }else {
        vast_rho.stRE[counter] <- paste0(rho.stRE.types[RhoConfig[3]+1], "-", rho.stRE.types[RhoConfig[4]+1])
      }
      
      
      #Setup File
      trial.dir <- paste0(working.dir,"/",n_x,"_bias.corr_",bias.correct)
      dir.create(trial.dir)
      trial.dir <- paste0(trial.dir, "/int_",vast_rho.int[counter], " stRE_",vast_rho.stRE[counter])
      dir.create(trial.dir)
    
    
      #=======================================================================
      ##### SNOWFALL CODE FOR PARALLEL #####
      sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
      sfExportAll() #Exportas all global variables to cores
      sfLibrary(TMB)  #Loads a package on all nodes
      sfLibrary(VAST)
      output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, RhoConfig=RhoConfig)
      sfStop()
    
      vast_est.output[[counter]] <- output
    }#next r
  }#next t
  
  #Dimensions for vast_est.output are 1) Trial knots, 2) Species
  # vast_est.output[[1:n.trial.knots]][[1:n.species]]
  
  #Create output directory
  dir.create(output.dir)
  save(vast_est.output, file=paste0(output.dir,"/vast_est.output.RData"))
  #Also save specifications
  vast_specs <- data.frame(vast_knots, vast_rho.int, vast_rho.stRE)
  write.csv(vast_specs, file=paste0(output.dir,"/vast_specs.csv"))
  
  #=======================================================================
  ##### DELETE UNNECESSARY FILE STRUCTURE #####
  #Must reset working directory
  # setwd(working.dir)
  # t <- 1
  # for(t in 1:n.trial.knots) {
  #   unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct), recursive=TRUE)
  # }#next t
  # 
  # time.2 <- date()
  # 
  # print(paste('### START:', time.1))
  # print(paste('### END:', time.2))
  
}else {
  load(paste0(output.dir,"/vast_est.output.RData"))
}














