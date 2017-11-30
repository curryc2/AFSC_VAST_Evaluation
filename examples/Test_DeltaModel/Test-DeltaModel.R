#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare VAST estimates with delta-GLMM
#
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 9.13.17
#
#Purpose: Evaluate whether differences between design-based and VAST (model-based) indices and uncertainty, is a function of the delta model structure. 
#
#
#==================================================================================================
#NOTES:

#==================================================================================================
#TIMING:
# [1] "### START: Wed Nov 29 15:35:28 2017"
# [1] "### END: Wed Nov 29 20:42:35 2017"



##==================================================================================================
require(parallel)
require(snowfall)
require(tidyverse)
require(ggthemes)
require(VAST)
require(TMB)

source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")
source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_DeltaModel")

#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

#Limit species included
species.list <- species.list[species.list$include=='Y',]
#Remove EBS_SHELF Arrowtooth
species.list <- species.list[species.list$survey!='EBS_SHELF',]

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
trial.knots <- c(100,500,1000)
n.trial.knots <- length(trial.knots)

#Trial RANDOM EFECTS SPECIFICATIONS specifications
trial.RE.names <- c('None','Spatial_Only','SpatioTemporal_Only','Full')
n.trial.RE <- length(trial.RE.names)
trial.RE <- matrix(c(0,0,0,0,
                     1,0,1,0,
                     0,1,0,1,
                     1,1,1,1), nrow=4, ncol=4, byrow=TRUE)

#Boolean for bias correction
bias.correct <- FALSE
#=======================================================================
##### Run VAST model  #####
Version <- "VAST_v2_8_0"
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
# FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#Output Directory Name
output.dir <- paste0(working.dir,"/output_bias.correct_",bias.correct)
dir.create(output.dir)


#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####

s <- 1
# for(s in 1:n.species) {
wrapper_fxn <- function(s, n_x, FieldConfig) {
  
  #Define file for analyses
  DateFile <- paste0(trial.dir,"/",species.list$survey[s],"_",species.list$name[s],"/")
  
  dir.create(DateFile)
  
  #Define species.codes
  species.codes <- species.list$species.code[s]
  survey <- species.list$survey[s]
  
  #=======================================================================
  ##### READ IN DATA AND BUILD VAST INPUT #####
  #  NOTE: this will create the DateFile
  VAST_input <- create_VAST_input(species.codes=species.codes, combineSpecies=FALSE,
                                  lat_lon.def=lat_lon.def, save.Record=FALSE,
                                  Method=Method, grid_size_km=grid_size_km, n_x=n_x,
                                  Kmeans_Config=Kmeans_Config,
                                  strata.limits=strata.limits, survey=survey,
                                  DateFile=DateFile,
                                  FieldConfig=FieldConfig, RhoConfig=RhoConfig,
                                  OverdispersionConfig=OverdispersionConfig,
                                  ObsModel=ObsModel, Options=Options, Version=Version)
  
  
  
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
  # Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
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
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List", "TmbList", "Obj")#, "Save")#, "Opt", "Report")
  
  #========================================================================
  setwd(home.dir)
  ##### RETURN SECTION #####
  out <- NULL
  out$vast_est <- vast_est
  out$Opt <- Opt
  return(out)
} 


#=======================================================================
##### Loop Through Trial Knots  #####
vast_est.output <- vector('list', length=(n.trial.knots * n.trial.RE))
vast_knots <- vector(length=(n.trial.knots * n.trial.RE))
vast_RE <- vector(length=(n.trial.knots * n.trial.RE))

if(do.estim==TRUE) {
  
  
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
    for(r in 1:n.trial.RE) {
      print(paste('#### Trial RE Scenario Number',r,'of',n.trial.RE))
      
      #Specify Spatial and Spatio-temporal RE
      FieldConfig <- trial.RE[r,]
      names(FieldConfig) <- c('Omega1','Epsilon1','Omega2','Epsilon2')
      
      #Record
      vast_RE[counter] <- trial.RE.names[r]
      vast_knots[counter] <- n_x
      
      #Setup File
      trial.dir <- paste0(working.dir,"/",n_x,"_bias.corr_",bias.correct)
      dir.create(trial.dir)
      trial.dir <- paste0(trial.dir, "/RE_",vast_RE[counter])
      dir.create(trial.dir)
      
      
      #=======================================================================
      ##### TEST WRAPPER FUNCTION #####
      # wrapper_fxn(s=1, n_x=n_x, FieldConfig=FieldConfig)
      
      #=======================================================================
      ##### SNOWFALL CODE FOR PARALLEL #####
      sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
      sfExportAll() #Exportas all global variables to cores
      sfLibrary(TMB)  #Loads a package on all nodes
      sfLibrary(VAST)
      output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, FieldConfig=FieldConfig)
      sfStop()

      vast_est.output[[counter]] <- output

      #For Update
      # output <- vast_est.output[[counter]]
      # save(output, file=paste0(output.dir, "/testVAST_output_",counter,".RData"), compression_level=9)
      saveRDS(output, file=paste0(output.dir, "/VAST_output_",counter,".rds"))
      #Real
      
      
      counter <- counter+1
    }#next r
  }#next t
  
  #Create output directory
  #Also save specifications
  vast_specs <- data.frame(vast_knots, vast_RE)
  write.csv(vast_specs, file=paste0(output.dir,"/vast_specs.csv"))
  
  #=======================================================================
  ##### DELETE UNNECESSARY FILE STRUCTURE #####
  #Must reset working directory
  setwd(working.dir)
  t <- 1
  for(t in 1:n.trial.knots) {
    unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct), recursive=TRUE)
  }#next t
  
  time.2 <- date()
  
  print(paste('### START:', time.1))
  print(paste('### END:', time.2))
  
}else {
  
  
  #Old
  # load(paste0(output.dir,"/vast_est.output.RData"))
  
  specs <- read.csv(paste0(output.dir,"/vast_specs.csv"), header=TRUE, stringsAsFactors=FALSE)
  n.specs <- nrow(specs)
  
  for(i in 1:n.specs) {
    print(i)
    vast_est.output[[i]] <- readRDS(file=paste0(output.dir, "/VAST_output_",i,".rds"))
    vast_knots[i] <- specs$vast_knots[i]
    vast_RE[i] <- specs$vast_RE[i]
  }#next i
}

