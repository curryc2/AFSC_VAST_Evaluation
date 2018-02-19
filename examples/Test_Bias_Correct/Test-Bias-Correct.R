#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare VAST estimates with and without BIAS CORRECTION
#
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 11.30.17
#
#Purpose: Determine if index scale is more sensitive sensitive to: a) spatial complexity (number of knots), or b) bias correction.
#
#
#==================================================================================================
#NOTES:
# 1)  Error in checkForRemoteErrors(val) : 8 nodes produced errors; first error: Memory allocation fail in function 'MakeADHessObject2' 
#  Suggests cannot be done in parallel, will attempt in series.

# 2) Running in series (32 GB Ram) stopped with bias.corr=TRUE, knots=400, s=4 (AI Walleye Pollock)
#==================================================================================================
#TIMING:




##==================================================================================================
require(parallel)
require(snowfall)
require(tidyverse)
require(ggthemes)
require(VAST)
require(TMB)
require(viridis)

source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")
source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Bias_Correct")

#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

#Limit species included
species.list <- species.list[species.list$include=='Y',]
#Remove EBS_SHELF Arrowtooth
species.list <- species.list[species.list$survey!='EBS_SHELF',]

#Remove AI Species
# species.list <- species.list[species.list$survey!='EBS_SHELF',]


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
<<<<<<< HEAD
trial.knots <- c(100,250)
=======
trial.knots <- seq(from=100, to=300, by=100)
>>>>>>> 0e8350bf30f3f8b4813d76d9588fbcd64cfafec7
n.trial.knots <- length(trial.knots)

#Trial RANDOM EFECTS SPECIFICATIONS specifications
trial.bias.correct <- c(FALSE,TRUE)
n.trial.bias.correct <- length(trial.bias.correct)


#Boolean for bias correction
# bias.correct <- FALSE
#=======================================================================
##### Run VAST model  #####
Version <- "VAST_v4_0_0"
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
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#Output Directory Name
output.dir <- paste0(working.dir,"/Testing")
dir.create(output.dir)


#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####

s <- 1
# for(s in 1:n.species) {
wrapper_fxn <- function(s, n_x, bias.correct) {
  
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
  cleanup_VAST_file(DateFile=DateFile, Version=Version) #No longer necessary as we are deleting everything at the end
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List", "TmbList", "Obj","Report")#, "Save")#, "Opt", "Report")
  
  dyn.unload("VAST_v2_8_0")
  # is.loaded()
  
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
vast_est.output <- vector('list', length=(n.species*n.trial.knots*n.trial.bias.correct))
vast_knots <- vector(length=(n.species*n.trial.knots*n.trial.bias.correct))
vast_bias.correct <- vector(length=(n.species*n.trial.knots*n.trial.bias.correct))
vast_species <- vector(length=(n.species*n.trial.knots*n.trial.bias.correct))

if(do.estim==TRUE) {
  
  
  time.1 <- date()
  
  #Counter for knots by rho
  counter <- 1
  s <- 1
  for(s in 1:n.species) {
    
    t <- 1
    for(t in 1:n.trial.knots) {
      print(paste('### Trial Species Number',s,'of',n.species))
      print(paste('## Trial Knot Number',t,'of',n.trial.knots))
      print(paste('# Trial Knots:',trial.knots[t]))
      #Specify trial observation model
    
      #Specify knots
      n_x <- trial.knots[t]
    
      r <- 1
      for(r in 1:n.trial.bias.correct) {
        print(paste('#### Trial Bias Correct',r,'of',n.trial.bias.correct))
      
<<<<<<< HEAD
        bias.correct <- trial.bias.correct[r]
        #Record
        vast_bias.correct[counter] <- bias.correct
        vast_knots[counter] <- n_x
        vast_species[counter] <- s
        #Setup File
        trial.dir <- paste0(working.dir,"/",n_x,"_bias.corr_",bias.correct)
        dir.create(trial.dir)
=======
      #=======================================================================
      ##### TEST WRAPPER FUNCTION #####
      output <- NULL
      output <- wrapper_fxn(s=s, n_x=n_x, bias.correct=bias.correct)
>>>>>>> 0e8350bf30f3f8b4813d76d9588fbcd64cfafec7
      
        #=======================================================================
        ##### TEST WRAPPER FUNCTION #####
        output <- wrapper_fxn(s=1, n_x=n_x, bias.correct=bias.correct)
      
        #=======================================================================
        ##### SNOWFALL CODE FOR PARALLEL #####
        # sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
        # sfExportAll() #Exportas all global variables to cores
        # sfLibrary(TMB)  #Loads a package on all nodes
        # sfLibrary(VAST)
        # output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, bias.correct=bias.correct)
        # sfStop()
        # 
        # vast_est.output[[counter]] <- output
      
        #For Update
        # output <- vast_est.output[[counter]]
        # save(output, file=paste0(output.dir, "/testVAST_output_",counter,".RData"), compression_level=9)
        saveRDS(output, file=paste0(output.dir, "/VAST_output_",counter,".rds"))
      
        counter <- counter+1
      }#next r
    }#next t
  }#next s  
  #Create output directory
  #Also save specifications
  vast_name <- paste0(species.list$survey[vast_species],"_",
                      species.list$name[vast_species])
  vast_specs <- data.frame(vast_species, vast_knots, vast_bias.correct,
                             vast_name)
  
  write.csv(vast_specs, file=paste0(output.dir,"/vast_specs.csv"))
  
  #=======================================================================
  ##### DELETE UNNECESSARY FILE STRUCTURE #####
  #Must reset working directory
  setwd(working.dir)
  t <- 1
  for(t in 1:n.trial.knots) {
    r <- 1
    for(r in 1:n.trial.bias.correct) {
      unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",trial.bias.correct[r]), recursive=TRUE)
    }#next r
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
    vast_bias.correct[i] <- specs$vast_bias.correct[i]
  }#next i
}

#=====================================================
# Gather Data
# vast.list <- NULL
# aic.list <- NULL
# aic.vect <- vector(length=0)
# converge.vect <- vector(length=0)
# 
# #Load dataset to determine which years to include
# goa.yrs <- sort(unique(load_RACE_data(species.codes=30420,
#                                       combineSpecies=FALSE, survey='GOA')$Year))
# ai.yrs <- sort(unique(load_RACE_data(species.codes=30420,
#                                      combineSpecies=FALSE, survey='AI')$Year))
# 
# 
# i <- 1
# for(i in 1:n.specs) {
#   s <- 1
#   for(s in 1:n.species) {
#     #Species Information
#     temp.species <- species.list$name[s]
#     temp.survey <- species.list$survey[s]
#     temp.name <- paste0(temp.survey,": ",temp.species)
#     
#     #Determine Survey years (Currently only GOA and AI)
#     if(temp.survey=='GOA') {
#       temp.yrs <- goa.yrs
#     }else {
#       temp.yrs <- ai.yrs
#     }
#     
#     #Get VAST model index
#     temp.list <- vast_est.output[[i]][[s]]$vast_est[c(1,4,6)]
#     
#     #Calculate CV
#     CV <- temp.list$SD_mt/temp.list$Estimate_metric_tons
#     
#     #Determine which are survey years
#     survey.year <- temp.list$Year %in% temp.yrs
#     
#     #Bind it
#     temp.list <- cbind(temp.list, CV, temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], vast_RE[i], survey.year)
#     
#     
#     #AIC and convergence
#     #Get AIC and convergence
#     temp.aic <- cbind(temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], vast_RE[i])
#     
#     aic.vect <- append(aic.vect, vast_est.output[[i]][[s]]$Opt$AIC)
#     converge.vect <- append(converge.vect, vast_est.output[[i]][[s]]$Opt$converge)
#     
#     #Combine to larger lists
#     vast.list <- rbind(vast.list, temp.list)
#     aic.list <- rbind(aic.list, temp.aic) 
#     
#     
#   }#Next species
#   
# }#next model configuration i
# 
# 
# #Add Design-based estimates
# db.list <- NULL
# s <- 1
# for(s in 1:n.species) {
#   #Species Information
#   temp.species <- species.list$name[s]
#   temp.survey <- species.list$survey[s]
#   temp.name <- paste0(temp.survey,": ",temp.species)
#   
#   #Determine Survey years (Currently only GOA and AI)
#   if(temp.survey=='GOA') {
#     temp.yrs <- goa.yrs
#   }else {
#     temp.yrs <- ai.yrs
#   }
#   
#   #Get design-based estimate
#   db_est <- calc_design_based_index(species.codes=species.list$species.code[s], survey=temp.survey)
#   
#   #Determine which are survey years
#   survey.year <- db_est$YEAR %in% temp.yrs  #All TRUE
#   
#   temp.name <- paste0(temp.survey,": ",temp.species)
#   temp.list <- cbind(db_est[,c(1,2,4,5)], temp.survey, temp.species, temp.name, "Design-based", 'Design-based','Design-based', survey.year)
#   #Add it
#   db.list <- rbind(db.list, temp.list)
# }#next s
# db.df <- data.frame(db.list)
# names(db.df) <- c('Year','Biomass','SD','CV','Survey','Species', 'Name','Model','Knots','RE','SurveyYear')
# 
# 
# #Add names
# vast.df <- data.frame(vast.list)
# names(vast.df) <- c('Year','Biomass','SD','CV','Survey','Species', 'Name','Model','Knots','RE','SurveyYear')
# 
# aic.df <- data.frame(aic.list, aic.vect, converge.vect)
# names(aic.df) <- c('Survey','Species','Name','Model','Knots','RE','AIC','Converge')
# aic.df$Converge <- as.factor(aic.df$Converge)
# 
# #Combine the lists
# survey.df <- rbind(vast.df,db.df)
# # survey.df$Knots <- factor(survey.df$Knots, ordered=TRUE)
# #=====================================================
# #PLOT IT OUT
# c('Walleye pollock','Pacific cod')
# -which(output.df$Species %in% c('Walleye pollock','Pacific cod',rockfish))
# 
# 
# ###### GOA Rockfish #####
# rockfish <- c('Pacific ocean perch','Northern rockfish','Harlequin rockfish')
# 
# plot.list <- survey.df[survey.df$Survey==survey &
#                          survey.df$Species %in% rockfish &
#                          survey.df$SurveyYear==TRUE,]
# # yrs.surv <- sort(unique(plot.list$Year[plot.list$Model=='Design-based']))
# # plot.list <- plot.list[plot.list$Year %in% yrs.surv,]
# 
# #Remove 2001 from design-based results because of incomplete sampling
# # if(survey=='GOA') { plot.list <- plot.list[-which(plot.list$Year==2001 & plot.list$Model=='Design-based'),] }
# plot.list$Biomass <- plot.list$Biomass/1e3
# 
# plot.list.v <- data.frame(plot.list[plot.list$Model!='Design-based',])
# plot.list.db <-  data.frame(plot.list[plot.list$Model=='Design-based', -which(names(plot.list)=='RE')])
# 
# g <- ggplot(plot.list[plot.list$Model!='Design-based',], 
#             aes(x=Year, y=Biomass, color=RE, lty=Model, ymin=0)) +
#   theme_dark() +
#   theme(legend.position='right') +
#   geom_line(lwd=1.5) +
#   geom_line(data=plot.list[plot.list$Model=='Design-based',], color='black') +
#   geom_point(data=plot.list[plot.list$Model=='Design-based',], show.legend=FALSE, colour='black') +
#   # geom_line(data=plot.list.db, color='black') +
#   # geom_point(data=plot.list.db, show.legend=FALSE, colour='black') +
#   # facet_wrap(~Species, scales='free', ncol=3) +
#   facet_grid(Species ~ Knots, scales='free') +
#   labs(list(y='Biomass (thousands of metric tonnes)')) +
#   ggtitle(paste(survey, 'Survey')) +
#   scale_color_viridis(discrete=TRUE) 
# 
# 
# g
# 
# 
