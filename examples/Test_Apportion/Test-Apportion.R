#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Apportionment Between VAST and
#                                                              RE Model for GOA
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 5.11.17
#
#Purpose: Evaluate apportionment (Eastern, Central, Western) for Gulf of Alaska RACE
#           Bottom Trawl indices of abundance. 
#             1) Evaluate across three alternative knot numbers
#             2) Several specifications for autocorrelation in EC and PCR intercept
#
#
#==================================================================================================
#NOTES:
#  a) Spiny Dogfish removed from evaluation because design based index is 0 for Western GOA in 2013
#  
#==================================================================================================
#TIMING:
##==================================================================================================

require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)
require(R2admb)
require(reshape2)
require(gridExtra)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")

source("R/run-RE-model.r")

home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Apportion")


#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

#Limit species included
species.list <- species.list[species.list$include=='Y',]
#Limit to GOA
species.list <- species.list[species.list$survey=='GOA',]
species.list <- species.list[species.list$name!='Spiny dogfish',]
n.species <- nrow(species.list)

#Create
species.series <- c(1:n.species)

#=======================================================================
##### CONTROL SECTION #####
#Specify process error terms for RE model
n_PE <- 3
PE_vec <- c(1:3)
#####

#Number of cores to use
n.cores <- detectCores()-1

#Boolean for running estimation models
do.estim <- FALSE

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
# trial.rho <- t(read.csv('Data/Test-Autoregressive-Input.csv', header=TRUE, stringsAsFactors=FALSE)[,-c(1:2)])
trial.rho <- matrix(c(2,2,0,0,
                      4,4,0,0),ncol=4, nrow=2, byrow=TRUE)
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
strata.limits <- data.frame(STRATA = c("Western","Central",'Eastern'),
                            west_border = c(-Inf, -159, -147),
                            east_border = c(-159, -147, Inf))


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
dir.create(output.dir)


#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####

s <- 9 #Spiny dogfish
# for(s in 1:n.species) {
wrapper_fxn <- function(s, n_x, RhoConfig, n_PE, PE_vec) {
  
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
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List",
     "TmbList", "Obj")#, "Save")#, "Opt", "Report")
  #========================================================================
  setwd(home.dir)
  #========================================================================
  ##### CALL RE MODEL #####
  #Calculate Separate Design-based indices for GOA
  temp.west <- calc_design_based_index(species.codes=species.codes, survey=survey, reg.area="WESTERN GOA")
  temp.cent <- calc_design_based_index(species.codes=species.codes, survey=survey, reg.area="CENTRAL GOA")  
  temp.east <- calc_design_based_index(species.codes=species.codes, survey=survey, reg.area="EASTERN GOA")
  
  #Bring Together (except 2001 when Eastern GOA was not surveyed)
  input.yrs <- sort(unique(temp.west$YEAR))
  #Remove 2001
  input.yrs <- input.yrs[-which(input.yrs==2001)]
  
  #Index values
  input.idx <- data.frame(temp.west$Biomass[temp.west$YEAR %in% input.yrs],
                          temp.cent$Biomass[temp.cent$YEAR %in% input.yrs],
                          temp.east$Biomass[temp.east$YEAR %in% input.yrs])
  names(input.idx) <- c("Western","Centeral","Eastern")
  input.idx <- as.matrix(input.idx)
  
  #Index CV
  input.cv <- data.frame(temp.west$CV[temp.west$YEAR %in% input.yrs],
                         temp.cent$CV[temp.cent$YEAR %in% input.yrs],
                         temp.east$CV[temp.east$YEAR %in% input.yrs])
  
  names(input.cv) <- c("Western","Centeral","Eastern")
  input.cv <- as.matrix(input.cv)
  
  #Copy, compile, call ADMB-RE model
  biomA <- run_RE_model(input.yrs, input.idx, input.cv, DateFile, home.dir, n_PE=n_PE, PE_vec=PE_vec)
    
  #========================================================================
  setwd(home.dir)
  
  ##### RETURN SECTION #####
  out <- NULL
  out$vast_est <- vast_est
  out$Opt <- Opt
  out$biomA <- biomA
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
      #Update Knot List
      vast_knots[counter] <- n_x
      
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
      output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, RhoConfig=RhoConfig, n_PE=n_PE, PE_vec=PE_vec)
      sfStop()
      
      #Add to list
      vast_est.output[[counter]] <- output
      #Save Object for storage
      saveRDS(output, file=paste0(output.dir,"/VAST_output_",counter,".rds"))
      #Update Counter
      counter <- counter + 1
    }#next r
  }#next t
  
  #Dimensions for vast_est.output are 1) Trial knots, 2) Species
  # vast_est.output[[1:n.trial.knots]][[1:n.species]]
  
  #Create output directory
  # dir.create(output.dir)
  # save(vast_est.output, file=paste0(output.dir,"/vast_est.output.RData"))
  #Also save specifications
  vast_specs <- data.frame(vast_knots, vast_rho.int, vast_rho.stRE)
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
  vast_est.output <- vector('list', length=(n.trial.knots * n.trial.rho))
  # vast_est.output <- array('list', dim=c(n.trial.knots,n.trial.rho))
  counter <- 1
  t <- 1
  for(t in 1:n.trial.knots) {
    r <- 1
    for(r in 1:n.trial.rho) {
      vast_est.output[[counter]] <- readRDS(file=paste0(output.dir,"/VAST_output_",counter,".rds"))
      counter <- counter + 1
    }#next r
  }#next t
  # load(paste0(output.dir,"/vast_est.output.RData"))
}

#Plot the output
#  First example is using not fully specified output


# plot.dat <- vast_est.output[[3]] #500 RW
# plot.dat <- vast_est.output[[4]] #500 AR
# plot.dat <- vast_est.output[[1]] #100 RW
#Determine years


#Loop through species
s <- 1
for(s in 1:n.species) {
  yrs <- sort(unique(plot.dat[[s]]$vast_est$Year))
  n.yrs <- length(yrs)

  indices <- c('Western','Central','Eastern')
  n.indices <- length(indices)

  temp.species <- species.list$name[s]

  #Random Effects Model
  temp.re <- plot.dat[[s]]$biomA
  prop.re <- temp.re
  y <- 1
  for(y in 1:n.yrs) {
    prop.re[y,] <- prop.re[y,]/sum(temp.re[y,])
  }
  dimnames(prop.re) <- list(yrs, indices)

  #VAST
  dat.vast <- plot.dat[[s]]$vast_est
  #Convert to a matrix
  temp.vast <- matrix(nrow=n.yrs, ncol=n.indices, dimnames=list(yrs, indices))
  i <- 1
  for(i in 1:n.indices) {
    temp.vast[,i] <- dat.vast$Estimate_metric_tons[dat.vast$Fleet==i]
  }
  prop.vast <- temp.vast
  y <- 1
  for(y in 1:n.yrs) {
    prop.vast[y,] <- prop.vast[y,]/sum(temp.vast[y,])
  }

  #Create the grand list
  list.re <- melt(prop.re)
  model <- 'ADMB-RE'
  list.re <- cbind(list.re,model,temp.species)

  list.vast <- melt(prop.vast)
  model <- 'VAST'
  list.vast <- cbind(list.vast,model,temp.species)
  #Combine
  if(s==1) {
    list.all <- rbind(list.re, list.vast)
  }else {
    list.all <- rbind(list.all, list.re, list.vast)
  }
}

names(list.all) <- c('Year','Region','value','Model','Species')



#===========================================================
#Plotting Rockfish
plot.rockfish <- FALSE
if(plot.rockfish==TRUE) {
#plot
rockfish <- c('Northern rockfish','Pacific ocean perch','Harlequin rockfish')
n.rockfish <- length(rockfish)
list.rf <- list.all[list.all$Species %in% rockfish,]



r <- 1
for(r in 1:n.rockfish) {

  g <- ggplot(list.rf[list.rf$Species==rockfish[r],], aes(x=Year, y=value, fill=Region)) +
                    theme_gray() +
                    geom_area(position='stack', alpha=0.75) +
                    facet_wrap(~Model, ncol=1) +
                    ggtitle(paste('Gulf of Alaska:',rockfish[r]))
  ggsave(paste0(output.dir,"/",rockfish[r],".png"), height=5, width=6, dpi=500, units='in')
}

#Together

g2 <- ggplot(list.rf, aes(x=Year, y=value, fill=Region)) +
        theme_gray() +
        geom_area(position='stack', alpha=0.75) +
        facet_grid(Model~Species)

ggsave(paste0(output.dir,"/Rockfish Apport.png"), plot=g2, height=5, width=8, dpi=500, units='in')

}


#===========================================================
#Plotting GOA Pollock
plot.GOA.pollock <- TRUE
if(plot.GOA.pollock==TRUE) {
  #plot
  specs <- c('Walleye pollock')
  temp.list <- list.all[list.all$Species %in% specs,]
  
  
  

    g <- ggplot(temp.list, aes(x=Year, y=value, fill=Region)) +
           theme_gray() +
           geom_area(position='stack', alpha=0.75) +
           facet_wrap(~Model, ncol=1) +
           ggtitle(paste('Gulf of Alaska:',specs[1]))
    ggsave(paste0(output.dir,"/",specs[1],"_100_RW.png"), height=5, width=6, dpi=500, units='in')
  # }
  
  #Together
  
  # g2 <- ggplot(list.rf, aes(x=Year, y=value, fill=Region)) +
  #   theme_gray() +
  #   geom_area(position='stack', alpha=0.75) +
  #   facet_grid(Model~Species)
  # 
  # ggsave(paste0(output.dir,"/Rockfish Apport.png"), plot=g2, height=5, width=8, dpi=500, units='in')
  
}

















