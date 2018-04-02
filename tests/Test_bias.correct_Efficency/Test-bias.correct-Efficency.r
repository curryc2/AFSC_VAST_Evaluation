#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Testing of Alternatives for Reducing RAM Load for Estimation with Epsilon-estimator for Bias Correction and Large (>300) Knots
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.1.2018
#
#Purpose: Example implementation of VAST model for GOA Pollock
#  bias.corr.option - controls which method is used.
#     Alternatives:  1 = Orignal implementation (error @ 300kt + )
#                    2 = Use nsplit to break up epsilon-estimator bias correction processes - slow
#                    3 = Manually estimte without bias correction then re-estimate the bias corrected index - faster
#                    4 = Use NEW option in TMBhelper::Optimize() to specify which variables are bias corrected.
#
#==================================================================================================
#NOTES:
#Memory ProfilingJust


# OPTION #2: Using nsplit=200
# [1] "bias.correct: TRUE"
# [1] "n_x: 1000"
# [1] "START: Mon Mar 05 10:18:23 2018"
# [1] "END: Mon Mar 05 17:28:07 2018"

# OPTION #3: Using Jim's new method - Restimate of index with bias corr. 
# [1] "bias.correct: TRUE"
# [1] "n_x: 1000"
# [1] "START: Tue Mar 06 10:02:39 2018"
# [1] "END: Tue Mar 06 11:24:57 2018"

# OPTION #4: New TMBHelper implementation - nsplit=200
# [1] "bias.correct: TRUE"
# [1] "n_x: 1000"
# [1] "START: Wed Mar 14 14:58:43 2018"
# [1] "END: Wed Mar 14 22:22:01 2018"
# If no nsplit specification, model has same memory allocation failure

#   Requires development version of TMB: devtools::install_github("kaskr/adcomp/TMB")
#     Checkup

#==================================================================================================
 source("R/create-VAST-input.r")
 source("R/create-Data-Geostat.r")
 source("R/load-RACE-data.r")
 source("R/plot-VAST-output.r")
 source("R/cleanup-VAST-file.r")
 source("R/run-RE-model.r") 

require(VAST)
require(TMB)
require(tidyverse)

#=======================================================================
##### SETUP INPUT DATA #####

#Generate a dataset
species.codes <- 21740
lat_lon.def <- "start"

survey <- "GOA"

bias.correct <- TRUE

#SPATIAL SETTINGS
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 25
n_x <- c(100, 250, 500, 1000, 2000)[4] # Number of stations
Kmeans_Config <- list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"))


#DERIVED OBJECTS
Version <-  "VAST_v4_0_0"

#SPECIFY METHOD FOR DEALING WITH BIAS.CORRECTION
bias.corr.option <- 4

###########################
trial.file <- paste0(getwd(),"/tests/Test_bias.correct_Efficency")

#MODEL SETTINGS
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal
# ObsModel = c(2, 0) #Gamma
# ObsModel = c(1, 1) #Poisson-Process Link function approximating Tweedie distribution

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)


DateFile <- paste0(trial.file,"/Option_",bias.corr.option,"/GOA Pollock knots_",n_x," bias.correct_", bias.correct, " Rho_",RhoConfig[1],RhoConfig[2],RhoConfig[3],RhoConfig[4],"/")


#=======================================================================
##### READ IN DATA AND BUILD vAST INPUT #####
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
Extrapolation_List <- VAST_input$Extrapolation_List #Becomes zeros for non-GOA
# head(Extrapolation_List$a_el)
# head(Extrapolation_List$Area_km2_x)
# head(Extrapolation_List$Data_Extrap)

#=======================================================================
##### RUN VAST #####



#Build TMB Object
#  Compilation may take some time
TmbList <- VAST::Build_TMB_Fn(TmbData = TmbData, RunDir = DateFile,
                                Version = Version, RhoConfig = RhoConfig, loc_x = Spatial_List$loc_x,
                                Method = Method)
Obj <- TmbList[["Obj"]]

start.time <- date()

#================================================
#TESTING OPTIMIZATION 1: Original Call
if(bias.corr.option==1) {
  Opt <- NULL
  Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                          upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                          bias.correct = bias.correct)
}

#================================================
#TESTING OPTIMIZATION 2: Updated call with nsplit to reduce memory load and 
#                        allow running bias.cor with kt > ~300
if(bias.corr.option==2) {
  nsplit <- 200
  Opt <- NULL
  Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                             upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                             bias.correct = bias.correct,
                             bias.correct.control=list(nsplit=nsplit, split=NULL, sd=FALSE))
}
#================================================
#TESTING OPTIMIZATION 3: New Alternative Following Jim's Suggestion
#  Should limit bias correction to single vector of interst: index
#  Conducts estimation without bias correction then re-estimates bias corrected index parameter.
if(bias.corr.option==3) {
  nsplit <- 200
  Opt <- NULL
  
  Opt = TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                            upper = TmbList[["Upper"]], getsd = TRUE,
                            savedir = DateFile, bias.correct=FALSE )

  # First SD run
  h <- optimHess(Opt$par, Obj$fn, Obj$gr)
  SD = sdreport( obj=Obj, par.fixed=Opt$par, hessian.fixed=h )

  # Determine indices
  BiasCorrNames = c("Index_cyl")
  Which = which( rownames(summary(SD,"report")) %in% BiasCorrNames )
  Which = split( Which, cut(seq_along(Which), nsplit) )
  Which = Which[sapply(Which,FUN=length)>0]


  # Repeat SD with indexing
  SD = sdreport( obj=Obj, par.fixed=Opt$par, hessian.fixed=h, bias.correct=TRUE, 
                 bias.correct.control=list(sd=FALSE, split=Which, nsplit=NULL) )
  
}
#================================================
#TESTING OPTIMIZATION 4: Updated TMB Implementation Where Transformed Variables Are Specified
if(bias.corr.option==4) {
  nsplit <- 100
  Opt <- NULL
  
  Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, 
                             savedir=DateFile, bias.correct=bias.correct, newtonsteps=1,  
                             # bias.correct.control=list(sd=FALSE, vars_to_correct="Index_cyl") )
  bias.correct.control=list(split=NULL, sd=FALSE, nsplit=nsplit, vars_to_correct="Index_cyl") )
   
}
#================================================


end.time <- date()
#Save output
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

#========================================================================
##### DIAGNOSTIC AND PREDICTION PLOTS #####
plot_VAST_output(Opt, Report, DateFile, survey, TmbData, Data_Geostat, Extrapolation_List, Spatial_List)

#========================================================================
##### CLEAN UP MODEL FILES #####
# cleanup_VAST_file(DateFile, Version=Version)

print(paste('bias.correct:',bias.correct))
print(paste('n_x:',n_x))
print(paste('nsplit:',nsplit))
print(paste('START:',start.time))
print(paste('END:',end.time))


#========================================================================
##### APPORTIONMENT #####














# 
# 
# #========================================================================
# ##### DIAGNOSTIC PLOTS #####
# 
# #Plot spatial distribution of data
# SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List = Extrapolation_List,
#                                       Spatial_List = Spatial_List, Data_Geostat = Data_Geostat,
#                                       PlotDir = DateFile)
# 
# #Diagnostics for Encounter Probability
# #  "Diag--Encounter_prob"
# Enc_prob = SpatialDeltaGLMM::Check_encounter_prob(Report = Report,
#                                                   Data_Geostat = Data_Geostat,
#                                                   DirName = DateFile)
# 
# #Diagnostics for positive-catch-rate component
# Q = SpatialDeltaGLMM::QQ_Fn(TmbData = TmbData, Report = Report,
#                             FileName_PP = paste0(DateFile, "Posterior_Predictive.jpg"),
#                             FileName_Phist = paste0(DateFile, "Posterior_Predictive-Histogram.jpg"),
#                             FileName_QQ = paste0(DateFile, "Q-Q_plot.jpg"),
#                             FileName_Qhist = paste0(DateFile, "Q-Q_hist.jpg"))
# 
# 
# #Diagnostics for plotting residuals on a map
# 
# 
# MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region,
#                                                    "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
#                                                    "Extrapolation_List"=Extrapolation_List )
# 
# #Which Years to Include
# Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
# Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
# 
# #Or just include years with observations
# 
# #Plot Pearson Residuals
# #  Look for spatial patterns-- indication of "overshrinking"
# #  Creates "maps--" files
# SpatialDeltaGLMM:::plot_residuals(Lat_i = Data_Geostat[,"Lat"], Lon_i = Data_Geostat[, "Lon"], TmbData = TmbData,
#                                   Report = Report, Q = Q, savedir = DateFile, MappingDetails = MapDetails_List[["MappingDetails"]],
#                                   PlotDF = MapDetails_List[["PlotDF"]], MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
#                                   Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]],
#                                   FileName = DateFile, Year_Set = Year_Set, Years2Include = Years2Include,
#                                   Rotate = MapDetails_List[["Rotate"]], Cex = MapDetails_List[["Cex"]],
#                                   Legend = MapDetails_List[["Legend"]], zone = MapDetails_List[["Zone"]],
#                                   mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8)
# 
# 
# 
# 
# #========================================================================
# ##### MODEL OUTPUT PLOTS #####
# 
# #Direction of "geometric anisotropy"
# SpatialDeltaGLMM::PlotAniso_Fn(FileName = paste0(DateFile,"Aniso.png"),
#                                Report = Report, TmbData = TmbData)
# 
# #Density Surface for Each Year -- "Dens"
# SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set = c(3),
#                                       MappingDetails = MapDetails_List[["MappingDetails"]],
#                                       Report = Report, Sdreport = Opt$SD, PlotDF = MapDetails_List[["PlotDF"]],
#                                       MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
#                                       Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]],
#                                       FileName = DateFile, Year_Set = Year_Set, Years2Include = Years2Include,
#                                       Rotate = MapDetails_List[["Rotate"]], Cex = MapDetails_List[["Cex"]],
#                                       Legend = MapDetails_List[["Legend"]], zone = MapDetails_List[["Zone"]],
#                                       mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8,
#                                       plot_legend_fig = TRUE)
# 
# 
# 
# #Generate Index of Abundance
# 
# Index = SpatialDeltaGLMM::PlotIndex_Fn(DirName = DateFile,
#                                        TmbData = TmbData, Sdreport = Opt[["SD"]],
#                                        Year_Set = Year_Set,
#                                        Years2Include = Years2Include, 
#                                        use_biascorr = TRUE)
# 
# idx <- Index$Table
# 
# 
# dev.off()
# #Plotting 
# 
# yrs.surv <- Year_Set[Years2Include]
# x.lim <- c(min(yrs.surv), max(yrs.surv))
# up.sd <- idx$Estimate_metric_tons + idx$SD_mt
# low.sd <- idx$Estimate_metric_tons - idx$SD_mt
# y.lim <- c(min(low.sd), max(up.sd))
# 
# loc.yrs <- which(idx$Year %in% yrs.surv)
# 
# 
# plot(x=NULL, y=NULL, xlim=x.lim, ylim=y.lim, ylab='Survey Estimate (metric Tons)', xlab='Year',
#      main='Gulf of Alaska\nNorthern Rockfish Survey Index')
# 
# polygon(x=c(yrs.surv, rev(yrs.surv)), y=c(low.sd[loc.yrs],rev(up.sd[loc.yrs])), col='lightblue', border=FALSE)
# lines(x=yrs.surv, y=idx$Estimate_metric_tons[loc.yrs], col='red')
# points(x=yrs.surv, y=idx$Estimate_metric_tons[loc.yrs], pch=21, bg='red')
# grid(col='black')
# 
# 
# #Center of gravity and range expansion/contraction
# #  For some reason I can't actually change the years to plot 
# SpatialDeltaGLMM::Plot_range_shifts(Report = Report,
#                                     TmbData = TmbData, Sdreport = Opt[["SD"]], Znames = colnames(TmbData$Z_xm),
#                                     PlotDir = DateFile, Year_Set = Year_Set)