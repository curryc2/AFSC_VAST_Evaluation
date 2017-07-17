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
require(ggthemes)


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

alpha <- 1

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
trial.knots <- c(100)
n.trial.knots <- length(trial.knots)

#Trial AUTOREGRESSIVE specifications
#Note starts at 0
# rho.int.types <- c('Fixed_Effect','Independent_Among_Years','Random_Walk','Constant_Intercept','Autoregressive')
rho.int.types <- c('FE','IaY','RW','CI','AR')
# rho.stRE.types <- c('Independent_Among_Years',NA,'Random_Walk',NA,'Autoregressive')
rho.stRE.types <- c('IaY',NA,'RW',NA,'AR')

#Read in Autoregressive Input
# trial.rho <- t(read.csv('Data/Test-Autoregressive-Input.csv', header=TRUE, stringsAsFactors=FALSE)[,-c(1:2)])
trial.rho <- matrix(c(1,1,0,0,
                      2,2,0,0,
                      4,4,0,0,
                      1,1,4,4,
                      2,2,4,4),ncol=4, nrow=5, byrow=TRUE)
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
Version <- "VAST_v2_5_0"
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

s <- 1 #Spiny dogfish
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
vast_est.output <- vector('list', length=(n.trial.knots * n.trial.rho))
vast_knots <- vector(length=(n.trial.knots * n.trial.rho))
vast_rho.int <- vector(length=(n.trial.knots * n.trial.rho))
vast_rho.stRE <- vector(length=(n.trial.knots * n.trial.rho))
  
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
  #Old
  # load(paste0(output.dir,"/vast_est.output.RData"))
  
  #New
  specs <- read.csv(paste0(output.dir,"/vast_specs.csv"), header=TRUE, stringsAsFactors=FALSE)
  n.specs <- nrow(specs)
  
  for(i in 1:n.specs) {
    print(i)
    vast_est.output[[i]] <- readRDS(file=paste0(output.dir, "/VAST_output_",i,".rds"))
    vast_knots[i] <- specs$vast_knots[i]
    vast_rho.int[i] <- specs$vast_rho.int[i]
    vast_rho.stRE[i] <- specs$vast_rho.stRE[i]
  }#next i
}

#Plot the output
#  First example is using not fully specified output


# plot.dat <- vast_est.output[[3]] #500 RW
# plot.dat <- vast_est.output[[4]] #500 AR
# plot.dat <- vast_est.output[[1]] #100 RW
#Determine years


#Loop through species
#
re.list <- NULL
vast.list <- NULL
aic.list <- NULL
aic.vect <- vector(length=0)
converge.vect <- vector(length=0)
maxGrad.vect <- vector(length=0)

#Determine true survey years for flag
goa.yrs <- sort(unique(load_RACE_data(species.codes=30420, survey='GOA')$Year))


s <- 1
for(s in 1:n.species) {
  #Species Information
  temp.species <- species.list$name[s]
  temp.survey <- species.list$survey[s]
  temp.name <- paste0(temp.survey,": ",temp.species)

  #VAST
  i <- 1
  for(i in 1:n.specs) { 
    
    #Update Convergence and AIC estimates
    aic.vect <- append(aic.vect, vast_est.output[[i]][[s]]$Opt$AIC)
    converge.vect <- append(converge.vect, vast_est.output[[i]][[s]]$Opt$converge)
    maxGrad.vect <- append(maxGrad.vect, max(abs(vast_est.output[[i]][[s]]$Opt$diagnostics$final_gradient)))
    
    yrs <- sort(unique(vast_est.output[[i]][[s]]$vast_est$Year))
    n.yrs <- length(yrs)
    
    survey.year <- vast_est.output[[i]][[s]]$vast_est$Year %in% goa.yrs
    
    indices <- c('Western','Central','Eastern')
    n.indices <- length(indices)

    #VAST
    dat.vast <- vast_est.output[[i]][[s]]$vast_est
    #Convert to a matrix
    temp.vast <- matrix(nrow=n.yrs, ncol=n.indices, dimnames=list(yrs, indices))
    j <- 1
    for(j in 1:n.indices) {
      temp.vast[,j] <- dat.vast$Estimate_metric_tons[dat.vast$Fleet==j]
    }
    prop.vast <- temp.vast
    y <- 1
    for(y in 1:n.yrs) {
      prop.vast[y,] <- prop.vast[y,]/sum(temp.vast[y,])
    }
    
    temp.RhoConfig <- paste(vast_rho.int[i], "+",vast_rho.stRE[i])
    
    #Create the grand list
    temp.vast <- melt(prop.vast)
    model <- 'VAST'
    temp.vast <- cbind(temp.vast, model, temp.survey, temp.species, temp.name, vast_knots[i], vast_rho.int[i], vast_rho.stRE[i], temp.RhoConfig, survey.year)

    #Skeleton in AIC/convergence data frame
    temp.aic <- cbind(temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], vast_rho.int[i], vast_rho.stRE[i], temp.RhoConfig)
    
    #Combine to larger lists
    #VAST list
    vast.list <- rbind(vast.list, temp.vast)
    #AIC
    aic.list <- rbind(aic.list, temp.aic) 
    
  }#next i
  
  #ADMB-RE
  #Random Effects Model
  temp.re <- vast_est.output[[i]][[s]]$biomA
  prop.re <- temp.re
  y <- 1
  for(y in 1:n.yrs) {
    prop.re[y,] <- prop.re[y,]/sum(temp.re[y,])
  }
  dimnames(prop.re) <- list(yrs, indices)
  
  #Create larger list
  temp.re <- melt(prop.re)
  model <- 'ADMB-RE'
  temp.re <- cbind(temp.re, model, temp.survey, temp.species, temp.name, FALSE, FALSE, FALSE, "ADMB-RE", survey.year)
  #Randome-effects list
  re.list <- rbind(re.list, temp.re)
  
}#next s

#Add names
re.df <- data.frame(re.list)
names(re.df) <- c('Year','Region','value','Model','Survey','Species', 'Name','Knots','Rho_Intercept','Rho_stRE','RhoConfig','SurveyYear')

vast.df <- data.frame(vast.list)
names(vast.df) <- c('Year','Region','value','Model','Survey','Species', 'Name','Knots','Rho_Intercept','Rho_stRE','RhoConfig','SurveyYear')

#Bind Together
output.df <- rbind(re.df, vast.df)

aic.df <- data.frame(aic.list, aic.vect, converge.vect, maxGrad.vect)
names(aic.df) <- c('Survey','Species','Name','Model','Knots','Rho_Intercept','Rho_stRE','RhoConfig','AIC','Converge','maxGradient')
aic.df$Converge <- as.factor(aic.df$Converge)

#===========================================================
#Remove the FE + IaY Model
#Unobserved estimates are unreliable
output.df <- output.df[output.df$RhoConfig!='FE + IaY',]
aic.df <- aic.df[aic.df$RhoConfig!='FE + IaY',]


#===========================================================
#Plot the Max Gradients
g <- ggplot(aic.df, aes(x=RhoConfig, y=maxGradient, color=Knots)) +
       theme_gray() +
       geom_point(alpha=0.5) +
       facet_wrap(~Species, scales='free') +
       theme(axis.text.x=element_text(angle=45, hjust=1))
g

#===========================================================
#Plotting Rockfish
    
rockfish <- c('Northern rockfish','Pacific ocean perch','Harlequin rockfish')
n.rockfish <- length(rockfish)
temp.df <- output.df[output.df$Species %in% rockfish,]
    
temp.knots <- 100
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,], 
              aes(x=Year, y=value, fill=Region)) +
       # theme_gray() +
       # theme_fivethirtyeight()+
       # theme_solarized()+
       # theme_economist() +
       theme_dark()+
       theme(legend.position='top') +
       # theme_tufte() +
       scale_fill_colorblind() +
       geom_area(position='stack', alpha=alpha) +
       facet_grid(Species~RhoConfig) +
       # facet_grid(RhoConfig~Species) +
       ggtitle(paste('Gulf of Alaska: Rockfish'), subtitle=paste('Knots:',temp.knots)) +
       ylab('Proportion')

g
ggsave(paste0(output.dir,'/Rockfish ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)

temp.knots <- 500
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,],
              aes(x=Year, y=value, fill=Region)) +
  theme_dark() +
  theme(legend.position='top') +
  scale_fill_colorblind() +
  ylab('Proportion') +
  geom_area(position='stack', alpha=alpha) +
  facet_grid(Species~RhoConfig) +
  ggtitle(paste('Gulf of Alaska: Rockfish'), subtitle=paste('Knots:',temp.knots))

g
ggsave(paste0(output.dir,'/Rockfish ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)

#===========================================================
#Plotting Pollock and Pcod
temp.df <- output.df[which(output.df$Species %in% c('Walleye pollock','Pacific cod')),]

temp.knots <- 100
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,],
              aes(x=Year, y=value, fill=Region)) +
  theme_dark() +
  theme(legend.position='top') +
  scale_fill_colorblind() +
  ylab('Proportion') +
  geom_area(position='stack', alpha=alpha) +
  facet_grid(Species~RhoConfig) +
  ggtitle(paste('Gulf of Alaska: Pollock and Cod'), subtitle=paste('Knots:',temp.knots))

g
ggsave(paste0(output.dir,'/Pollock Cod ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)


temp.knots <- 500
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,], aes(x=Year, y=value, fill=Region)) +
  theme_dark() +
  theme(legend.position='top') +
  scale_fill_colorblind() +
  ylab('Proportion') +
  geom_area(position='stack', alpha=alpha) +
  facet_grid(Species~RhoConfig) +
  ggtitle(paste('Gulf of Alaska: Pollock and Cod'), subtitle=paste('Knots:',temp.knots))

g
ggsave(paste0(output.dir,'/Pollock Cod ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)

#===========================================================
#Plotting Pollock and Pcod

temp.df <- output.df[-which(output.df$Species %in% c('Walleye pollock','Pacific cod',rockfish)),]

temp.knots <- 100
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,],
            aes(x=Year, y=value, fill=Region)) +
  theme_dark() +
  theme(legend.position='top') +
  scale_fill_colorblind() +
  ylab('Proportion') +
  geom_area(position='stack', alpha=alpha) +
  facet_grid(Species~RhoConfig) +
  ggtitle(paste('Gulf of Alaska: Other Species'), subtitle=paste('Knots:',temp.knots))

g
ggsave(paste0(output.dir,'/Other Species ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)


temp.knots <- 500
g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,], aes(x=Year, y=value, fill=Region)) +
  theme_dark() +
  theme(legend.position='top') +
  scale_fill_colorblind() +
  ylab('Proportion') +
  geom_area(position='stack', alpha=alpha) +
  facet_grid(Species~RhoConfig) +
  ggtitle(paste('Gulf of Alaska: Other Species'), subtitle=paste('Knots:',temp.knots))

g
ggsave(paste0(output.dir,'/Other Species ', temp.knots,'kt.png'), plot=g, height=8, width=8, units='in', dpi=500)


# cbind(vast_est.output[[1]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[2]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[3]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[4]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[5]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[6]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[7]][[7]]$vast_est$Estimate_metric_tons,
#       vast_est.output[[8]][[7]]$vast_est$Estimate_metric_tons)


