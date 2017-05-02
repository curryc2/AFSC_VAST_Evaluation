#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Generating design-based estimates for comparison
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 3.30.17
#
#Purpose: To explore efficient estimation of design-based index for RACE bottom trawl data.
#
#
#
#==================================================================================================
#NOTES:
#
#  a) Should be implemented as a function which calls RACE data read fxn.
#  b) Should return point estimate and SD by year.
#  c) Selection of strata should be clear
#==================================================================================================
require(VAST)
require(TMB)

source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")
#=======================================================================
##### CONTROL SECTION #####


working.dir <- getwd()


#NEW: Run all species... to lazy for parallel right now.


#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)
#Limit to those included
species.list <- species.list[species.list$include=='Y',]

n.species <- nrow(species.list)
species.series <- c(1:n.species)



#=======================================================================
##### Run VAST model  #####
lat_lon.def <- "start" 

#SPATIAL SETTINGS
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
n_x = c(100, 250, 500, 1000, 2000)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


#SET SRATIFICATOIN
#Basic - Single Area
strata.limits <- data.frame(STRATA = c("All_areas"))


#DERIVED OBJECTS
Version <-  "VAST_v2_4_0"

#Bias correction
bias.correct <- FALSE

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

start.time <- date()

s <- 14
# for(s in 1:n.species) {
  print(paste(s, 'of', n.species))
  species.codes <- species.list$species.code[s]
  survey <- species.list$survey[s]
  spec.name <- species.list$name[s]
  
  #=======================================================================
  ##### Calculate design-based estimate  #####
  db_est <- calc_design_based_index(species.codes=species.codes, survey=survey)
  
  ###########################
  DateFile=paste0(getwd(),'/examples/Test_DesignBased_Estimator/')
  #Create input
  VAST_input <- create_VAST_input(species.codes=species.codes, lat_lon.def=lat_lon.def, save.Record=TRUE,
                                  Method=Method, grid_size_km=grid_size_km, n_X=n_X,
                                  Kmeans_Config=Kmeans_Config,
                                  strata.limits=NULL, survey=survey,
                                  DateFile=DateFile,
                                  FieldConfig, RhoConfig, OverdispersionConfig,
                                  ObsModel, Options)

  #Unpack
  TmbData <- VAST_input$TmbData
  Data_Geostat <- VAST_input$Data_Geostat
  Spatial_List <- VAST_input$Spatial_List
  Extrapolation_List <- VAST_input$Extrapolation_List


  #Build TMB Object
  #  Compilation may take some time
  TmbList <- VAST::Build_TMB_Fn(TmbData = TmbData, RunDir = DateFile,
                                Version = "VAST_v2_4_0", RhoConfig = RhoConfig, loc_x = Spatial_List$loc_x,
                                Method = Method)
  Obj <- TmbList[["Obj"]]


  Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                             upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                             bias.correct = bias.correct)

  #=======================================================================
  ##### Plot comparison of design-based and model-based estimates #####

  Year_Set <- seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include <- which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  years <- Year_Set[Years2Include]

  vast_est <- get_VAST_index(TmbData=TmbData, Sdreport=Opt[["SD"]], bias.correct=bias.correct, Data_Geostat=Data_Geostat)
  #Limit to years with observations
  vast_est <- vast_est[Years2Include,]
  
  pdf(paste0(DateFile, survey, "_", spec.name, "_", TmbData$n_x, " knots.pdf"), height=8, width=8)
  
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,1))
  
  
  y.lim <- c(0, max(vast_est$Estimate_metric_tons+2*vast_est$SD_mt,
                    db_est$Biomass+2*db_est$SD))
  x.lim <- c(min(Year_Set), max(Year_Set))

  #Plot it out
  plot(x=NULL, y=NULL, xlim=x.lim, ylim=y.lim, xlab='Year', ylab='Estimator',
         main=paste0(spec.name, ' n=', TmbData$n_x, ' knots'))
  grid(col='black')

  legend('bottomleft', legend=c('Design-based', 'VAST'), fill=c('blue', 'red'), ncol=2, bg='white')

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
  
  #Cleanup unnecessary VAST Program Files after minimization
  if(s==n.species) { cleanup_VAST_file(DateFile=DateFile, Version=Version) }
  
  #Need to reset working director PlotIndex_Fn seems to change to DateFile
  setwd(working.dir)

} #end wrapper function

# dev.off()


end.time <- date()

print(paste('START:',start.time))
print(paste('END:',end.time))
