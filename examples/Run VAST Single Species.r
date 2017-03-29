#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Example Script
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.28.17
#
#Purpose: Example of How to Use helper functions
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
source('R/create-VAST-input.r')

require(VAST)
require(TMB)

#=======================================================================
##### SETUP INPUT DATA #####

#Generate a dataset
species.codes=c(30420) #Rockfish
lat_lon.def="mean"

#SPATIAL SETTINGS
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
n_x = c(100, 250, 500, 1000, 2000)[1] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"),
                            west_border = c(-Inf),
                            east_border = c(Inf))


#DERIVED OBJECTS
Region = "Gulf_of_Alaska"
###########################
DateFile=paste0(getwd(),'/examples/VAST_output/')

#MODEL SETTINGS
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Delta Model
# ObsModel = c(1, 1) #Poisson-Process Link function approximating Tweedie distribution

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#=======================================================================
##### READ IN DATA AND BUILD vAST INPUT #####




out <- create_VAST_input(species.codes=c(30420), lat_lon.def="mean", save.Record=TRUE,
                                     Method="Mesh", grid_size_km=25, n_X=250,
                                     Kmeans_Config=list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 ),
                                     strata.limits=NULL, Region="Gulf_of_Alaska",
                                     DateFile=DateFile,
                                     FieldConfig, RhoConfig, OverdispersionConfig,
                                     ObsModel, Options)



#Unpack
TmbData <- out$TmbData
Data_Geostat <- out$Data_Geostat
Spatial_List <- out$Spatial_List



#=======================================================================
##### RUN VAST #####



#Build TMB Object
#  Compilation may take some time
TmbList <- VAST::Build_TMB_Fn(TmbData = TmbData, RunDir = DateFile,
                                Version = "VAST_v2_4_0", RhoConfig = RhoConfig, loc_x = Spatial_List$loc_x,
                                Method = Method)
Obj <- TmbList[["Obj"]]


Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                          upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                          bias.correct = FALSE)
#Save output
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))


