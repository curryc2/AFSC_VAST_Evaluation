#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Testing Observation Model Specification
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.30.17
#
#Purpose: To test impact of changing distribution on observation model for positive catch rates,
#           across species.
#  1) Loop through obs model specifications for positive catch rate distribution.
#  2) Run vAST model in parallel across species. 
#  3) Initially run with 100 knots.
#
#
#==================================================================================================
#NOTES:
#  a) Originally included Normal distribution for positive catch rates,
#       but produced unknown errors for some species.
#  b) Obs models that returned errors for some species:
#       [0] Normal
#       [5] Negative binomial
#       [7] Poisson
#  c) Conway-Maxwell Poisson seems to run without error, but is EXTREMELY SLOW
#  d) QQ_Fn() from SpatialDeltaGLMM was VERY SLOW.
#       Below code is implemented to generate the QQplot more efficiently, but 
#         currently only for gamma and lognormal obs models. 
#
#
#==================================================================================================
require(snowfall)
require(parallel)
require(ggplot2)
require(TMB)
require(TMBhelper)
require(VAST)
require(reshape)
require(foreach)

#Source necessary files
source("R/create-VAST-input.r")
source("R/cleanup-VAST-file.r")

#Create working directory
working.dir <- paste0(getwd(),'/examples/Test_Obs_Model')

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
do.estim <- FALSE

#Setup observation model trials
# trial.obs <- c(0:2,5:7)
trial.obs <- c(0:2,5,7)
n.trial.obs <- length(trial.obs)

#Names for trial Obs models
# names.trial.obs <- c('Normal','Lognormal','Gamma','Negative_binomial','Conway-Maxwell_Poisson','Poisson')
names.trial.obs <- c('Normal','Lognormal','Gamma','Negative_binomial','Poisson')
#=======================================================================
##### VAST MODEL SPECIFICATIONS #####

Version <- "VAST_v2_4_0"

lat_lon.def <- "mean"

#SPATIAL SETTINGS
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 25
n_x <- c(100, 250, 500, 1000, 2000)[1] # Number of stations
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
# for(s in 1:n.species) {
species_wrapper_fxn <- function(s) {
  
  #Define file for analyses
  DateFile <- paste0(trial.dir,"/",species.list$name[s],"/")
  
  #Define species.codes
  species.codes <- species.list$species.code[s]
  
  #=======================================================================
  ##### READ IN DATA AND BUILD VAST INPUT #####
  #  NOTE: this will create the DateFile
  
  VAST_input <- create_VAST_input(species.codes=species.codes, lat_lon.def=lat_lon.def, save.Record=TRUE,
                                  Method=Method, grid_size_km=grid_size_km, n_X=n_x,
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
  
  rm("VAST_input", "TmbData", "Data_Geosta", "Spatial_List", "Extrapolation_List",
     "TmbList", "Obj", "Opt", "Report", "Save")
  
  #========================================================================
  ##### RETURN SECTION #####
  # return(Opt$AIC)
  
} 

if(do.estim==TRUE) {
  

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
    # sfRemove(Save)
    # sfRemover(VAST_input)
    sfStop()
  
    #Output
    AICs[t,] <- unlist(rbind(output))
  
  }# next t
}

#=======================================================================
##### SUMMARY AND PLOTTING SECTION #####
#NOTE: Only Gamma and Lognormal worked for all stocks.

names.species <- species.list$name

#Compare AIC across all potential models
temp.obs <- c('Normal','Lognormal','Gamma','Negative_binomial','Poisson')
n.temp.obs <- length(temp.obs)

#Loop through and determine if model ran... and if so grab AIC
ran.mat <- array(dim=c(n.temp.obs, n.species), dimnames=list(temp.obs,names.species)) #Boolean matrix for whether model ran
aic.mat <- array(dim=c(n.temp.obs, n.species), dimnames=list(temp.obs,names.species))

i <- 1
for(i in 1:n.temp.obs) {
  j <- 1
  for(j in 1:n.species) {
    obs <- temp.obs[i]
    spec <- names.species[j]
    
    temp.dir <- paste0(working.dir,"/",obs,"/",spec)
    
    #Determine if Model ran
    if(file.exists(paste0(temp.dir,"/parameter_estimates.txt"))) {
      ran.mat[i,j] <- TRUE
      load(paste0(temp.dir,"/parameter_estimates.RData"))
      aic.mat[i,j] <- parameter_estimates$AIC
    }else {
      ran.mat[i,j] <- FALSE
      aic.mat[i,j] <- NA
    }
  }#next j
}#next i

#Plot as AIC bars
ran.list <- melt(ran.mat)
names(ran.list) <- c('Obs','Species','value')
aic.list <- melt(aic.mat)
names(aic.list) <- c('Obs','Species','value')

g <- ggplot(aic.list, aes(x=Species, y=value, fill=Obs)) +
       theme_gray() +
       geom_col(position='dodge') +
       theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       labs(list(y='AIC', fill='Obs. Error\nDist.'))
# g
ggsave(paste0(working.dir,"/figs/AIC Compare at 100 knots.png"), plot=g, height=7, width=7)

#COMPARE LOGNORMAL VS GAMMA
temp.obs <- c('Lognormal','Gamma')
n.temp.obs <- length(temp.obs)

#Create directory for qqplot output\

qq.dir <- paste0(working.dir,"/figs/QQplot")
dir.create(qq.dir)

pdf(paste0(working.dir, "/figs/QQ Plots.pdf"), height=8, width=6)

pow = function(a,b) a^b

par(mfrow=c(n.species/2,n.temp.obs), mar=c(2,2,2,0), oma=c(3,3,0,0))

j <- 1
for(j in 1:n.species) {
  i <- 1
  for(i in 1:n.temp.obs) {
    obs <- temp.obs[i]
    spec <- names.species[j]
    temp.dir <- paste0(working.dir,"/",obs,"/",spec)

    #Load Results
    load(paste0(temp.dir,"/Save.RData"))

    #Plot
    TmbData <- Save$TmbData
    Report <- Save$Report


     #Jim's original function. Time consuming...  so lets speed it up
     # Q <-  SpatialDeltaGLMM::QQ_Fn(TmbData = TmbData, Report = Report,
     #                               FileName_QQ = paste0(qq.dir, "/", spec, "_", obs, ".jpg"))
    
    
    #Taken from Jim's QQ_Fn code, simplified by Curry for lognorm/gamma and
    #  reduced processing time with vectorized format and foreach call. 
  
    Which = which(TmbData$b_i>0)
    Q = rep(NA, length(Which) ) # vector to track quantiles for each observation
    y = array(NA, dim=c(length(Which),1000))
    pred_y = var_y = rep(NA, length(Which) ) # vector to track quantiles for each observation

    #  Simulate quantiles for different distributions
    if(TmbData$ObsModel[1]==1){
      pred_y <- TmbData$a_i[Which] *exp(Report$P2_i[Which])
      
      # set.seed(1)
      # for(ObsI in 1:length(Which)){
      #   # pred_y[ObsI] = TmbData$a_i[Which[ObsI]] * exp(Report$P2_i[Which[ObsI]])
      #   y[ObsI,] = rlnorm(n=ncol(y), meanlog=log(pred_y[ObsI])-pow(Report$SigmaM[1],2)/2, sdlog=Report$SigmaM[1])   # Plotting in log-space
      #   # Q[ObsI] = plnorm(q=TmbData$b_i[Which[ObsI]], meanlog=log(pred_y[ObsI])-pow(Report$SigmaM[1],2)/2, sdlog=Report$SigmaM[1])
      # }
      # set.seed(1)
      #Run Parallel
      y <- foreach(ObsI=c(1:length(Which)), .combine='rbind') %do% {
        rlnorm(n=ncol(y), meanlog=log(pred_y[ObsI])-pow(Report$SigmaM[1],2)/2, sdlog=Report$SigmaM[1])
      }
  
      #Alternative in vector format
      Q <- plnorm(q=TmbData$b_i[Which], meanlog=log(pred_y)-pow(Report$SigmaM[1],2)/2, sdlog=Report$SigmaM[1])
      
      
    }
    if(TmbData$ObsModel[1]==2){
      pred_y <- TmbData$a_i[Which] * exp(Report$P2_i[Which])
      b <- pow(Report$SigmaM[1],2) * pred_y;
      #Loop it
      # set.seed(123)
      # for(ObsI in 1:length(Which)){
      #   y[ObsI,] = rgamma(n=ncol(y), shape=1/pow(Report$SigmaM[1],2), scale=b)
      # }
      # set.seed(123)
      #Foreach it
      t.y <- foreach(ObsI=c(1:length(Which)), .combine='rbind') %do% {
        rgamma(n=ncol(y), shape=1/pow(Report$SigmaM[1],2), scale=b)
      }
      
      Q <- pgamma(q=TmbData$b_i[Which], shape=1/pow(Report$SigmaM[1],2), scale=b)
    }
    #Plot it out
    Qtemp = na.omit(Q)
    Order = order(Qtemp)
    plot(x=seq(0,1,length=length(Order)), y=Qtemp[Order], main=paste0(obs, ": ", spec), xlab="", ylab="", type="l", lwd=3)
    abline(a=0,b=1)

  }#next i
  mtext('Empirical', side=2, outer=TRUE, font=2)
  mtext('Uniform', side=1, outer=TRUE, font=2)
}#next j

dev.off()









