#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Temporal Linkages Between
#                                                              Intercepts and Spatio-temporal RE
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 4.10.17
#
#Purpose: To explore VAST model-based index sensitivity to autoregressive specifications on:
#           1) Intercept
#           2) Spatio-temporal random effect
#           3) Fix the number of knots @ 500, to begin with and consider changing to c(100,500) for sensitivity
#           4) Given 
#
#
#
#==================================================================================================
#NOTES:
#  a) int_RW-FE stRE_IaY - Fails on EBS Arrowtooth flounder. 
#  b) As a result we have removed EBS_SHELF Arrowtooth from this example.
#  
#==================================================================================================
#TIMING:
#
#==================================================================================================
require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)
require(ggthemes)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Autoregressive")

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
trial.knots <- c(100,500)
n.trial.knots <- length(trial.knots)

#Trial AUTOREGRESSIVE specifications
#Note starts at 0
# rho.int.types <- c('Fixed_Effect','Independent_Among_Years','Random_Walk','Constant_Intercept','Autoregressive')
rho.int.types <- c('FE','IaY','RW','CI','AR')
# rho.stRE.types <- c('Independent_Among_Years',NA,'Random_Walk',NA,'Autoregressive')
rho.stRE.types <- c('IaY',NA,'RW',NA,'AR')

#Read in Autoregressive Input
trial.rho <- t(read.csv('Data/Test-Autoregressive-Input.csv', header=TRUE, stringsAsFactors=FALSE)[,-c(1:2)])
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
strata.limits <- data.frame(STRATA = c("All_areas"))

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

s <- 1
# for(s in 1:n.species) {
wrapper_fxn <- function(s, n_x, RhoConfig) {
  
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
  
  rm("VAST_input", "TmbData", "Data_Geostat", "Spatial_List", "Extrapolation_List",
     "TmbList", "Obj")#, "Save")#, "Opt", "Report")
  
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
      output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, RhoConfig=RhoConfig)
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
  
  #Dimensions for vast_est.output are 1) Trial knots, 2) Species
  # vast_est.output[[1:n.trial.knots]][[1:n.species]]
  
  #Create output directory
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


#=====================================================
# Gather Data
vast.list <- NULL
aic.list <- NULL
aic.vect <- vector(length=0)
converge.vect <- vector(length=0)

#Load dataset to determine which years to include
goa.yrs <- sort(unique(load_RACE_data(species.codes=30420, survey='GOA')$Year))
ai.yrs <- sort(unique(load_RACE_data(species.codes=30420, survey='AI')$Year))


i <- 1
for(i in 1:n.specs) {
  s <- 1
  for(s in 1:n.species) {
    #Species Information
    temp.species <- species.list$name[s]
    temp.survey <- species.list$survey[s]
    temp.name <- paste0(temp.survey,": ",temp.species)
    
    #Determine Survey years (Currently only GOA and AI)
    if(temp.survey=='GOA') {
      temp.yrs <- goa.yrs
    }else {
      temp.yrs <- ai.yrs
    }
    
    #Get VAST model index
    temp.list <- vast_est.output[[i]][[s]]$vast_est[c(1,4,6)]
    
    #Calculate CV
    CV <- temp.list$SD_mt/temp.list$Estimate_metric_tons
    
    #Determine which are survey years
    survey.year <- temp.list$Year %in% temp.yrs
    
    #Bind it
    temp.list <- cbind(temp.list, CV, temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], vast_rho.int[i], vast_rho.stRE[i], survey.year)
    
    
    #AIC and convergence
    #Get AIC and convergence
    temp.aic <- cbind(temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], vast_rho.int[i], vast_rho.stRE[i])

    aic.vect <- append(aic.vect, vast_est.output[[i]][[s]]$Opt$AIC)
    converge.vect <- append(converge.vect, vast_est.output[[i]][[s]]$Opt$converge)
    
    #Combine to larger lists
    vast.list <- rbind(vast.list, temp.list)
    aic.list <- rbind(aic.list, temp.aic) 
    
    
  }#Next species
  
}#next model configuration i



#Add names
vast.df <- data.frame(vast.list)
names(vast.df) <- c('Year','Biomass','SD','CV','Survey','Species', 'Name','Model','Knots','Rho_Intercept','Rho_stRE','SurveyYear')

head(vast.df)

aic.df <- data.frame(aic.list, aic.vect, converge.vect)
names(aic.df) <- c('Survey','Species','Name','Model','Knots','Rho_Intercept','Rho_stRE','AIC','Converge')
aic.df$Converge <- as.factor(aic.df$Converge)


#=====================================================
#PLOT IT OUT
g <- ggplot(aic.df, aes(x=Rho_Intercept, y=AIC, colour=Knots, shape=Converge)) +
       theme_gray() +
       geom_point(size=5, alpha=0.5) +
       facet_wrap(Survey~Species, scales='free')
g

# ggsave(paste0(output.dir,'/AIC Compare.png'), plot=g, height=7, width=9, units='in', dpi=500)

g <- ggplot(aic.df[aic.df$Survey=='GOA',], aes(x=Rho_Intercept, y=AIC, colour=Knots)) +
       theme_gray() +
       geom_point(size=5, alpha=0.5) +
       facet_wrap(~Species, scales='free') +
       ggtitle('Gulf of Alaska Survey') +
       theme(axis.text.x=element_text(angle=90, hjust=1))
g
ggsave(paste0(output.dir,'/AIC Compare_GOA.png'), plot=g, height=7, width=9, units='in', dpi=500)

g <- ggplot(aic.df[aic.df$Survey=='AI',], aes(x=Rho_Intercept, y=AIC, colour=Knots)) +
       theme_gray() +
       geom_point(size=5, alpha=0.5) +
       facet_wrap(~Species, scales='free') +
       ggtitle('Aleutian Islands Survey') +
       theme(axis.text.x=element_text(angle=90, hjust=1))
g
ggsave(paste0(output.dir,'/AIC Compare_AI.png'), plot=g, height=7, width=9, units='in', dpi=500)


#===================================
#Plotting Index Values

g <- ggplot(aic.df, aes(x=Knots, y=AIC, colour=Rho_Intercept)) +
  theme_gray() +
  # geom_point(size=5, alpha=0.5) +
  facet_wrap(Survey~Species, scales='free') +
  geom_jitter()
# g

g <- ggplot(aic.df, aes(x=Rho_Intercept, y=AIC, colour=Knots)) +
  theme_gray() +
  # geom_point(size=5, alpha=0.5) +
  geom_jitter()
  facet_grid(Species~Survey, scales='free')
# g


#Plot some trends

#Determine which years to include

vast.df$Knots <- as.factor(vast.df$Knots)

#species x knots
g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$SurveyYear==TRUE,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  geom_line() +
  geom_point() +
  facet_grid(Species~Knots, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')
  
g

#Break it down by species groups
goa.species <- as.vector(unique(vast.df$Species[vast.df$Survey=='GOA' & vast.df$SurveyYear==TRUE]))

goa.rockfish <- c('Pacific ocean perch', 'Northern rockfish', 'Harlequin rockfish')
goa.others <- goa.species[-which(goa.species %in% goa.rockfish)]

#Plot rockfish
g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$SurveyYear==TRUE & vast.df$Species%in%goa.rockfish,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~Knots, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g

ggsave(paste0(output.dir,'/Rockfish Trend Compare_GOA.png'), plot=g, height=9, width=8, units='in', dpi=500)
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/Rockfish Trend Compare_GOA.png"), plot=g, height=7, width=10, units='in', dpi=1e3)


#Plot others
g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$SurveyYear==TRUE & vast.df$Species%in%goa.others,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  # scale_color_colorblind()+
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # scale_color_solarized() +
  facet_grid(Species~Knots, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')
  # theme(legend.position='bottom')

g

ggsave(paste0(output.dir,'/Others Trend Compare_GOA.png'), plot=g, height=9, width=8, units='in', dpi=500)
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/Others Trend Compare_GOA.png"), plot=g, height=7.75, width=11, units='in', dpi=1e3)

#================
g <- ggplot(vast.df[vast.df$Survey=='AI' & vast.df$SurveyYear==TRUE,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  facet_grid(Species~Knots, scales='free') +
  ggtitle('Aleutian Islands Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g
ggsave(paste0(output.dir,'/Trend Compare_AI.png'), plot=g, height=9, width=8, units='in', dpi=500)
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/Trend Compare_AI.png"), plot=g, height=7, width=10, units='in', dpi=1e3)


#=================================================
g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$Knots==100 &
            vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
            vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/Survey_non Compare Rockfish_GOA.png"), plot=g, height=7, width=10, units='in', dpi=1e3)

g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.others & vast.df$Knots==100 &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/Survey_non Compare others_GOA.png"), plot=g, height=7.75, width=11, units='in', dpi=1e3)
#=================================================




g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$Knots==100,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g






g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$Knots==100,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g


g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$Knots==500 &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g

g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$SurveyYear==TRUE,],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Knots, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_wrap(Species~Rho_Intercept, scales='free', nrow=length(goa.rockfish)) +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g


#=====================================================================================================
g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.rockfish & vast.df$Knots==100,],# &
                      # vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      # vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g

g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Species%in%goa.others & vast.df$Knots==100 &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g


g <- ggplot(vast.df[vast.df$Survey=='AI' & vast.df$Knots==100 &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~SurveyYear, scales='free') +
  ggtitle('Aleutian Islands Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g

g <- ggplot(vast.df[vast.df$Survey=='AI' &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  facet_grid(Species~Knots, scales='free') +
  ggtitle('Aleutian Islands Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g


g <- ggplot(vast.df[vast.df$Survey=='GOA' & vast.df$Knots==100 &
                      vast.df$Rho_Intercept!='AR-FE' & vast.df$Rho_Intercept!='RW-FE' &
                      vast.df$Rho_Intercept!='FE-AR' & vast.df$Rho_Intercept!='FE-RW',],
            aes(x=Year, y=Biomass/1e3, group=Rho_Intercept, colour=Rho_Intercept, ymin=0)) +
  theme_gray() +
  # scale_color_colorblind() +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  # theme_fivethirtyeight() +
  # scale_color_colorblind() +
  # facet_grid(Species~Knots, scales='free') +
  facet_wrap(~Species, scales='free') +
  ggtitle('Gulf of Alaska Survey') +
  ylab('Biomass (thousands of metric tonnes)')

g



# g <- ggplot(vast.list[vast.list$Survey=='GOA',],
#             aes(x=Knots, y=Biomass, group=Rho_Intercept, fill=Rho_Intercept, ymin=0)) +
#   theme_gray() +
#   geom_boxplot() +
#   facet_wrap(~Species, scales='free')
# g

# g <- ggplot(vast.list[vast.list$Survey=='GOA' & vast.list$Rho_Intercept!='RW-FE',],
#             aes(x=Year, y=Biomass/1e3,  colour=Rho_Intercept, ymin=0)) +
#   theme_gray() +
#   geom_line() +
#   facet_grid(Species~Knots, scales='free')
# g

