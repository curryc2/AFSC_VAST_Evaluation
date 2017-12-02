#' Wrapper function generate VAST diagnostic and prediction plots
#'
#' @param Opt list output from Template Model Builder (TMB) fitting of VAST model
#' @param Report list of TMB specifications for VAST model
#' @param DateFile path for directory housing VAST model and output figures and objects
#' @param survey string indicating the survey for which data are being extracted: GOA, AI, EBS_SHELF, EBS_SLOPE
#' @param TmbData list of TMB inputs for VAST model
#' @param Data_Geostat dataframe containing RACE bottom trawl survey data
#' @param Extrapolation_List list object containing definition of for which density extrapolation is conducted when calculating indices
#' @param Spatial_List tagged list with spatial information needed for Data_Fn
#'
#' @return A series of output figures exploring model diagnostics, density predictions, and indices.
#' @export
plot_VAST_output <- function(Opt, Report, DateFile, survey, TmbData, Data_Geostat, Extrapolation_List, Spatial_List) {
  
  #Determine Correct area to allign with RACE survey
  if(survey %in% c("GOA","AI","EBS_SHELF",'EBS_SLOPE')) { 
    if(survey=="GOA") { Region <- "Gulf_of_Alaska"; area <- "GOA" }
    if(survey=="AI") { Region <- "Aleutian_Islands"; area <- "AI" }
    if(survey=="EBS_SHELF" | survey=="EBS_SLOPE") { Region <- "Eastern_Bering_Sea"; area <- "BS" }
    
  }else {
    stop(paste("survey is:",survey,", should be one of: GOA, AI, EBS_SHELF, EBS_SLOPE"))
  }
  
  #========================================================================
  ##### DIAGNOSTIC PLOTS #####
  
  #Plot spatial distribution of data
  SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List = Extrapolation_List,
                                        Spatial_List = Spatial_List, Data_Geostat = Data_Geostat,
                                        PlotDir = DateFile)
  
  #Diagnostics for Encounter Probability
  #  "Diag--Encounter_prob"
  Enc_prob = SpatialDeltaGLMM::Check_encounter_prob(Report = Report,
                                                    Data_Geostat = Data_Geostat,
                                                    DirName = DateFile)
  
  #Diagnostics for positive-catch-rate component
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
  # Year_Set = y
  # Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  
  #Plot Pearson Residuals
  #  Look for spatial patterns-- indication of "overshrinking"
  #  Creates "maps--" files
  SpatialDeltaGLMM::plot_residuals(Lat_i = Data_Geostat[,"Lat"], Lon_i = Data_Geostat[, "Lon"], TmbData = TmbData,
                                    Report = Report, Q = Q, savedir = DateFile, MappingDetails = MapDetails_List[["MappingDetails"]],
                                    PlotDF = MapDetails_List[["PlotDF"]], MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
                                    Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]],
                                    FileName = DateFile, Year_Set = Year_Set, Years2Include = Years2Include,
                                    Rotate = MapDetails_List[["Rotate"]], Cex = MapDetails_List[["Cex"]],
                                    Legend = MapDetails_List[["Legend"]], zone = MapDetails_List[["Zone"]],
                                    mar = c(0, 0, 2, 0), oma = c(3.5, 3.5, 0, 0), cex = 1.8)
  
  
  
  
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
  
  # idx <- Index$Table
  # 
  # 
  # #Plotting 
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
  # polygon(x=c(yrs.surv, rev(yrs.surv)), y=c(low.sd[loc.yrs],rev(up.sd[loc.yrs])), col='lightgray', border=FALSE)
  # lines(x=yrs.surv, y=idx$Estimate_metric_tons[loc.yrs], col='red')
  # points(x=yrs.surv, y=idx$Estimate_metric_tons[loc.yrs], pch=21, bg='red')
  # grid(col='black')
  
  
  #Center of gravity and range expansion/contraction
  #  For some reason I can't actually change the years to plot 
  SpatialDeltaGLMM::Plot_range_shifts(Report = Report,
                                      TmbData = TmbData, Sdreport = Opt[["SD"]], Znames = colnames(TmbData$Z_xm),
                                      PlotDir = DateFile, Year_Set = Year_Set)
  
  
}
