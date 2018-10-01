#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Gulf of Alaska Northern Rockfish for Comparison to spatialDeltaGLMM()
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 5.22.17
#
#Purpose: Example implementation of VAST model for GOA Northern rockfish, 
#           to be proposed for 2018 stock assessment
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")
source("R/plot-VAST-output.r")
source("R/cleanup-VAST-file.r")
source("R/run-RE-model.r") 

require(VAST)
require(TMB)

#=======================================================================
##### SETUP INPUT DATA #####

#Generate a dataset
species.codes <- c(30420)
combineSpecies <- FALSE

lat_lon.def <- "start"

survey <- "GOA"
Region <- 'Gulf_of_Alaska'

bias.correct <- TRUE

#SPATIAL SETTINGS
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 25
n_x <- 500 #c(100, 250, 500, 1000, 2000)[1] # Number of stations
Kmeans_Config <- list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"))


#DERIVED OBJECTS
Version <-  "VAST_v4_0_0"
###########################
trial.file <- paste0(getwd(),"/examples/Species_Specific_Case_Studies/Test_GOA_Northern_rockfish/")

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


DateFile <- paste0(trial.file,"GOA Northern rockfish knots_",n_x," bias.correct_", bias.correct, " Rho_",RhoConfig[1],RhoConfig[2],RhoConfig[3],RhoConfig[4],"/")
#=======================================================================
##### READ IN DATA AND BUILD vAST INPUT #####


#Extract Dusy data for comparison with Dana's version 10.12.17
temp.data <- load_RACE_data(species.codes=species.codes, combineSpecies=combineSpecies, survey=survey, writeCSV=FALSE)
write.csv(temp.data, file=paste0(trial.file,"/Northern rockfish Input Data Curry.csv"))


VAST_input <- create_VAST_input(species.codes=species.codes, combineSpecies=combineSpecies,
                                     lat_lon.def=lat_lon.def, save.Record=TRUE,
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


# Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
#                           upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
#                           bias.correct = bias.correct,
#                           bias.correct.control=list(sd=FALSE, nsplit=200, split=NULL,
#                                vars_to_correct="Index_cyl"),
#                           newtonsteps = 2)
#                           # bias.correct.control=list(nsplit=200, split=NULL, sd=FALSE))

if(bias.correct==FALSE) {
  Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                             upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                             bias.correct = bias.correct, newtonsteps=2)
}else {
  #NEW: Only Bias Correct Index
  Opt <- TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], 
                             upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, 
                             bias.correct=bias.correct, newtonsteps=2,
                             bias.correct.control=list(sd=TRUE, nsplit=NULL, split=NULL,
                                                       vars_to_correct="Index_cyl"))
}
#Save output
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

#========================================================================
##### DIAGNOSTIC AND PREDICTION PLOTS #####
# plot_VAST_output(Opt, Report, DateFile, survey, TmbData, Data_Geostat, Extrapolation_List, Spatial_List)

#========================================================================
##### CLEAN UP MODEL FILES #####
# cleanup_VAST_file(DateFile, Version=Version)


# #========================================================================
# ##### DIAGNOSTIC PLOTS #####
# 
#Plot spatial distribution of data
SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List = Extrapolation_List,
                                      Spatial_List = Spatial_List, Data_Geostat = Data_Geostat,
                                      PlotDir = DateFile)

#Diagnostics for Encounter Probability
#  "Diag--Encounter_prob"
Enc_prob = SpatialDeltaGLMM::Check_encounter_prob(Report = Report,
                                                  Data_Geostat = Data_Geostat,
                                                  DirName = DateFile)

#Diagnostics for positive-catch-rate component - WARNINGS
Q = SpatialDeltaGLMM::QQ_Fn(TmbData = TmbData, Report = Report,
                            FileName_PP = paste0(DateFile, "Posterior_Predictive.jpg"),
                            FileName_Phist = paste0(DateFile, "Posterior_Predictive-Histogram.jpg"),
                            FileName_QQ = paste0(DateFile, "Q-Q_plot.jpg"),
                            FileName_Qhist = paste0(DateFile, "Q-Q_hist.jpg"))


#Diagnostics for plotting residuals on a map


MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region,
                                                   "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
                                                   "Extrapolation_List"=Extrapolation_List )

#Which Years to Include
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

#Or just include years with observations

#Plot Pearson Residuals - NOT WORKING
#  Look for spatial patterns-- indication of "overshrinking"
#  Creates "maps--" files
# SpatialDeltaGLMM:::plot_residuals(Lat_i = Data_Geostat[,"Lat"], Lon_i = Data_Geostat[, "Lon"], TmbData = TmbData,
#                                   Report = Report, Q = Q, savedir = DateFile, MappingDetails = MapDetails_List[["MappingDetails"]],
#                                   PlotDF = MapDetails_List[["PlotDF"]], MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
#                                   Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]],
#                                   FileName = DateFile, Year_Set = Year_Set, Years2Include = Years2Include,
#                                   Rotate = MapDetails_List[["Rotate"]], Cex = MapDetails_List[["Cex"]],
#                                   Legend = MapDetails_List[["Legend"]], zone = MapDetails_List[["Zone"]],
#                                   mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8)




SpatialDeltaGLMM:::plot_residuals(Lat_i = Data_Geostat[,"Lat"], Lon_i = Data_Geostat[, "Lon"], TmbData = TmbData, 
                                  Report = Report, Q = Q, savedir = DateFile, MappingDetails = MapDetails_List[["MappingDetails"]], 
                                  PlotDF = MapDetails_List[["PlotDF"]], MapSizeRatio = MapDetails_List[["MapSizeRatio"]], 
                                  Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]], 
                                  FileName = DateFile, Year_Set = Year_Set, Rotate = 1,#MapDetails_List[["Rotate"]], 
                                  Cex = MapDetails_List[["Cex"]], Legend = MapDetails_List[["Legend"]], 
                                  zone = MapDetails_List[["Zone"]], mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8)

#========================================================================
##### MODEL OUTPUT PLOTS #####

#Direction of "geometric anisotropy"
SpatialDeltaGLMM::PlotAniso_Fn(FileName = paste0(DateFile,"Aniso.png"),
                               Report = Report, TmbData = TmbData)

#Density Surface for Each Year -- "Dens"
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set = c(3),
                                      MappingDetails = MapDetails_List[["MappingDetails"]],
                                      Report = Report, Sdreport = Opt$SD, PlotDF = MapDetails_List[["PlotDF"]],
                                      MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
                                      Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]],
                                      FileName = DateFile, Year_Set = Year_Set, Years2Include = Years2Include,
                                      Rotate = MapDetails_List[["Rotate"]], Cex = MapDetails_List[["Cex"]],
                                      Legend = MapDetails_List[["Legend"]], zone = MapDetails_List[["Zone"]],
                                      mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8,
                                      plot_legend_fig = TRUE)



#Generate Index of Abundance

Index = SpatialDeltaGLMM::PlotIndex_Fn(DirName = DateFile,
                                       TmbData = TmbData, Sdreport = Opt[["SD"]],
                                       Year_Set = Year_Set,
                                       Years2Include = Years2Include,
                                       use_biascorr = TRUE)

idx <- Index$Table


dev.off()
#Plotting

# Extract VAST Index =======================================================
source("R/get-VAST-index.r")
vast_est = get_VAST_index(TmbData=TmbData, Sdreport=Opt[["SD"]], bias.correct=bias.correct, Data_Geostat=Data_Geostat)
#Limit to years with observations
vast_est = vast_est[Years2Include,]

# Compare VAST and Design-based ============================================
#Calculate Design-Based
source("R/load-RACE-data.R")
source("R/calc-design-based-index.R")
db_est = calc_design_based_index(species.codes=species.codes, survey=survey)

# Plot 

png(paste0(DateFile,"/VAST DB Index Comparison.png"), height=8, width=9, units='in', res=500)  
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,1))

y.lim = c(0, max(vast_est$Estimate_metric_tons+2*vast_est$SD_mt,
                 db_est$Biomass+2*db_est$SD, na.rm=TRUE))
x.lim = c(min(Year_Set), max(Year_Set))

#Survey years to plot
years = Year_Set[Years2Include]

#Plot it out
plot(x=NULL, y=NULL, xlim=x.lim, ylim=y.lim, xlab='Year', ylab='Abundance (metric tonnes)',
     main=paste0(' n=', TmbData$n_x, ' knots'))
grid(col='black')

legend('top', legend=c('Design-based', 'VAST'), fill=c('blue', 'red'), ncol=2, bg='white')

#Design-based
polygon(x=c(years, rev(years)), y=c(db_est$Biomass+2*db_est$SD, rev(db_est$Biomass-2*db_est$SD)),
        border=FALSE, col=rgb(0,0,1, alpha=0.25))

polygon(x=c(years, rev(years)), y=c(db_est$Biomass+1*db_est$SD, rev(db_est$Biomass-1*db_est$SD)),
        border=FALSE, col=rgb(0,0,1, alpha=0.25))

lines(x=years, y=db_est$Biomass, lwd=2, col='blue')
points(x=years, y=db_est$Biomass, pch=21, bg='blue')

#VAST
polygon(x=c(years, rev(years)), y=c(vast_est$Estimate_metric_tons+2*vast_est$SD_mt,
                                    rev(vast_est$Estimate_metric_tons-2*vast_est$SD_mt)),
        border=FALSE, col=rgb(1,0,0, alpha=0.25))

polygon(x=c(years, rev(years)), y=c(vast_est$Estimate_metric_tons+1*vast_est$SD_mt, 
                                    rev(vast_est$Estimate_metric_tons-1*vast_est$SD_mt)),
        border=FALSE, col=rgb(1,0,0, alpha=0.25))

lines(x=years, y=vast_est$Estimate_metric_tons, lwd=2, col='red')
points(x=years, y=vast_est$Estimate_metric_tons, pch=21, bg='red')
dev.off()

