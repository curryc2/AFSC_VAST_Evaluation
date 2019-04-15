#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Apportionment Between VAST and
#                                                              RE Model for GOA - Running in Series to Permit Bias Corr
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
#  a) Spiny Dogfish removed from evaluation because design based index is 0 for Western GOA in 2013.
#  b) Big Skate Causes Problems... not diagnosed.
#  c) Non-convergence with intercepts estimted as IaY i.e. ==1, independent of spatio-temporal RE specs.
#  d) Must estimate spatio-temporal random effects, fails to converge without: flag.spatial <- c(TRUE) ONLY
#  e) Convergence failure when trying to do partial bias correction results in the 3-column error message. 
#  f) rho  4,4,4,4, fails on Northern Rockfish, Epsilon_rho1 (autocorr coeff for EP) has HIGH gradient.
#  g) GOA Harlequin Rockfish causes problem with AR - hessian not positive definite.
#      Epsilon_rho1, Epsilon_rho2, and Beta_rho2 all at upper limits.
#  h) ObsModel = c(1, 1) #Lognormal - Poisson-link Delta, 100 knots - Spiny Dogfish - rho=c(2,2,2,2)
#       ADMB Fail "Error: Invalid index 4 used for array range [1, 3] in "dvector& dmatrix::operator[] (int i)". 
#                     matrix bound exceeded -- row index too high"

#  i) Big Skate has trouble with AR intercepts, mixed convergence max grad 0.02
#==================================================================================================
#TIMING:
# 100 knots obs model c(1,1) - lognormal with poisson-link delta - no bias correction, spatial==TRUE only, and rho (2,2,2,2 and 4,4,4,4)
# [1] "### START: Tue Jan 29 09:38:47 2019"
# [1] "### END: Tue Jan 29 13:25:35 2019"

# [1] "### START: Sun Apr 14 22:50:32 2019"
# [1] "### END: Mon Apr 15 06:46:08 2019"

#Lognormal with Po
##==================================================================================================
# require(SpatialDeltaGLMM)
require(VAST)
require(TMB)
require(TMBhelper)
require(parallel)
require(snowfall)
require(ggplot2)
require(R2admb)
require(reshape2)
require(gridExtra)
require(ggthemes)
require(cowplot)
require(RANN)
require(PBSmodelling)
require(PBSmapping)
require(XML)
require(FishStatsUtils)


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
# species.list <- species.list[species.list$name!='Spiny dogfish' & species.list$name!='Big skate',]
species.list <- species.list[species.list$name!='Spiny dogfish',]

n.species <- nrow(species.list)

#Create
species.series <- c(1:n.species)

#=======================================================================
##### CONTROL SECTION #####
#Specify process error terms for RE model
# n_PE <- 3
# PE_vec <- c(1:3)

n_PE <- 1
PE_vec <- rep(1,3)
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
# trial.rho <- matrix(c(1,1,0,0,
#                       2,2,0,0,
#                       4,4,0,0,
#                       1,1,4,4,
#                       1,1,2,2,
#                       2,2,4,4),ncol=4, nrow=6, byrow=TRUE)

# trial.rho <- matrix(c(2,2,0,0,
#                       4,4,0,0,
#                       2,2,4,4,
#                       4,4,2,2),ncol=4, nrow=4, byrow=TRUE)

# trial.rho <- matrix(c(2,2,0,0,
#                       4,4,0,0,
#                       2,2,2,2,
#                       4,4,4,4),ncol=4, nrow=4, byrow=TRUE)

# trial.rho <- matrix(c(2,2,2,2,
#                       4,4,4,4
#                       ),ncol=4, nrow=2, byrow=TRUE)

# trial.rho <- matrix(c(2,2,2,2,
#                       4,4,2,2
#                       ),ncol=4, nrow=2, byrow=TRUE)

trial.rho <- matrix(c(0,0,2,2,
                      2,2,2,2,
                      4,4,2,2
                      ),ncol=4, nrow=3, byrow=TRUE)

n.trial.rho <- nrow(trial.rho)

flag.spatial <- c(TRUE)#,FALSE)
n.flag.spatial <- length(flag.spatial)
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
Version <- "VAST_v4_0_0"
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
# FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
# RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 1) #Lognormal - Poisson-link Delta
# ObsModel = c(1, 0) #Lognormal - Delta

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 0, Calculate_evenness = 0, Calculate_effective_area = 0,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#Output Directory Name
output.dir <- paste0(working.dir,"/output_bias.correct_",bias.correct, " n_PE-",n_PE)
dir.create(output.dir)


#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####

# s <- 1 #Spiny dogfish
# for(s in 1:n.species) {

#=======================================================================
##### Loop Through Trial Knots  ##### 
vast_est.output <- vector('list', length=(n.trial.knots * n.trial.rho * n.flag.spatial))
vast_knots <- vector(length=(n.trial.knots * n.trial.rho * n.flag.spatial))
vast_rho.int <- vector(length=(n.trial.knots * n.trial.rho * n.flag.spatial))
vast_rho.stRE <- vector(length=(n.trial.knots * n.trial.rho * n.flag.spatial))
vast_do.spatial <-  vector(length=(n.trial.knots * n.trial.rho * n.flag.spatial))

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
    # r <- 3
    for(r in 1:n.trial.rho) {
      #Specify intercepts and spatio-temporal variation across time
      RhoConfig <- trial.rho[r,]
      names(RhoConfig) <- c('Beta1','Beta2','Epsilon1','Epsilon2')
      
      f <- 1
      for(f in 1:n.flag.spatial) {
        do.spatial <- flag.spatial[f]
        # print(paste('f',f))
        #Turn ON/OFF spatial components
        if(do.spatial==TRUE) {
          #ON
          FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
        }else {
          #OFF
          FieldConfig = c(Omega1 = 0, Epsilon1 = 1, Omega2 = 0, Epsilon2 = 1)
        }
        vast_do.spatial[counter] <- do.spatial
        
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
        trial.dir <- paste0(trial.dir, "/int_",vast_rho.int[counter], 
                            " stRE_",vast_rho.stRE[counter], 
                            " do.spatial_",vast_do.spatial[counter])
        dir.create(trial.dir)
        
        
        #=======================================================================
        ##### SNOWFALL CODE FOR PARALLEL #####
        # sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
        # sfExportAll() #Exportas all global variables to cores
        # sfLibrary(TMB)  #Loads a package on all nodes
        # sfLibrary(VAST)
        # output <- sfLapply(species.series, fun=wrapper_fxn, n_x=n_x, RhoConfig=RhoConfig,
        #                    n_PE=n_PE, PE_vec=PE_vec, FieldConfig=FieldConfig)
        # #
        # # temp.out <- wrapper_fxn(s=1, n_x=n_x, RhoConfig=RhoConfig,
        # #                         n_PE=n_PE, PE_vec=PE_vec, FieldConfig=FieldConfig, bias.correct=bias.correct, species.list=species.list)
        # 
        # sfStop()
        
        output <- vector('list', length=n.species)
        s <- 1
        for(s in species.series) {
          # output[[s]]
          # wrapper_fxn <- function(ss, n_x, RhoConfig, n_PE, PE_vec, FieldConfig) {
            # require(TMB)
            # require(TMBhelper)
            # require(VAST)
            
            #Define file for analyses
            DateFile <- paste0(trial.dir,"/",species.list$survey[s],"_",species.list$name[s],"/")
            
            dir.create(DateFile, recursive=TRUE)
            
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
            
            # str(VAST_input)
            #Unpack
            TmbData <- VAST_input$TmbData
            Data_Geostat <- VAST_input$Data_Geostat
            Spatial_List <- VAST_input$Spatial_List
            Extrapolation_List <- VAST_input$Extrapolation_List
            
            # print(TmbData)
            
            print('Before create TmbList')
            #=======================================================================
            ##### RUN VAST #####
            #Build TMB Object
            #  Compilation may take some time - ERROR HERE
            TmbList <- VAST::Build_TMB_Fn(TmbData = TmbData, RunDir = DateFile,
                                          Version = Version, 
                                          Q_Config=FALSE, CovConfig=FALSE,
                                          RhoConfig = RhoConfig, loc_x = Spatial_List$loc_x,
                                          Method = Method)
            print('After create TmbList')
            
            Obj <- TmbList[["Obj"]]
            
            if(bias.correct==FALSE) {
              Opt <- TMBhelper::Optimize(obj = Obj, lower = TmbList[["Lower"]],
                                         upper = TmbList[["Upper"]], getsd = TRUE, savedir = DateFile,
                                         bias.correct = bias.correct, newtonsteps=1)
              # summary(Opt)
            }else {
              # #NEW: Only Bias Correct Index
              # Opt <- TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], 
              #                            upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, 
              #                            bias.correct=bias.correct, newtonsteps=1,
              #                            bias.correct.control=list(sd=TRUE, nsplit=10, split=NULL,
              #                                                      vars_to_correct="Index_cyl"))
              Opt <- TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], 
                                  upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, 
                                  bias.correct=bias.correct, newtonsteps=1)
            }
            #Save output
            # Report = Obj$report()
            # Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
            # save(Save, file=paste0(DateFile,"Save.RData"))
            
            #HESSIAN NOT POSITIVE DEFINITE
            # bs <- sdreport( obj=Obj, bias.correct=FALSE )  #
            # Opt$SD <- bs
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
               "TmbList", "Obj","Save","Report")#, "Save")#, "Opt", "Report")
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
          #   return(out)
          # } 
          
          #RETURN SECTION ============
          output[[s]] <- out
          
        }#next s
        
        #Add to list
        vast_est.output[[counter]] <- output
        #Save Object for storage
        saveRDS(output, file=paste0(output.dir,"/VAST_output_",counter,".rds"))
        #Update Counter
        # print(counter)
        counter <- counter + 1
      }#next f
      
    }#next r
  }#next t
  
  #Dimensions for vast_est.output are 1) Trial knots, 2) Species
  # vast_est.output[[1:n.trial.knots]][[1:n.species]]
  
  #Create output directory
  # dir.create(output.dir)
  # save(vast_est.output, file=paste0(output.dir,"/vast_est.output.RData"))
  #Also save specifications
  vast_specs <- data.frame(vast_knots, vast_rho.int, vast_rho.stRE, vast_do.spatial)
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
  
  #Need to re-set working directory
  setwd(home.dir)
  
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
    vast_do.spatial[i] <- specs$vast_do.spatial[i]
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
    converge.vect <- append(converge.vect, ifelse(vast_est.output[[i]][[s]]$Opt$Convergence_check=="There is no evidence that the model is not converged", TRUE, FALSE))
    # maxGrad.vect <- append(maxGrad.vect, max(abs(vast_est.output[[i]][[s]]$Opt$diagnostics$final_gradient)))  #Equivalent
    maxGrad.vect <- append(maxGrad.vect, vast_est.output[[i]][[s]]$Opt$max_gradient)
    
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
    temp.vast <- cbind(temp.vast, model, temp.survey, temp.species, temp.name, vast_knots[i], 
                       vast_rho.int[i], vast_rho.stRE[i], temp.RhoConfig, survey.year,
                       vast_do.spatial[i])
    
    #Skeleton in AIC/convergence data frame
    temp.aic <- cbind(temp.survey, temp.species, temp.name, 'VAST', vast_knots[i], 
                      vast_rho.int[i], vast_rho.stRE[i], temp.RhoConfig,
                      vast_do.spatial[i])
    
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
  temp.re <- cbind(temp.re, model, temp.survey, temp.species, temp.name, FALSE, FALSE, FALSE, "ADMB-RE", 
                   survey.year, NA)
  #Randome-effects list
  re.list <- rbind(re.list, temp.re)
  
}#next s

#Add names
re.df <- data.frame(re.list)
names(re.df) <- c('Year','Region','value','Model','Survey','Species', 'Name','Knots',
                  'Rho_Intercept','Rho_stRE','RhoConfig','SurveyYear','Est_Spatial_RE')

vast.df <- data.frame(vast.list)
names(vast.df) <- c('Year','Region','value','Model','Survey','Species', 'Name','Knots',
                    'Rho_Intercept','Rho_stRE','RhoConfig','SurveyYear','Est_Spatial_RE')

#Bind Together
output.df <- rbind(re.df, vast.df)

aic.df <- data.frame(aic.list, aic.vect, converge.vect, maxGrad.vect)

names(aic.df) <- c('Survey','Species','Name','Model','Knots',
                   'Rho_Intercept','Rho_stRE','RhoConfig','Est_Spatial_RE',
                   'AIC','Converge','maxGradient')
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
  geom_point() +
  # facet_grid(Est_Spatial_RE~Species, scales='free') +
  facet_wrap(~Species, ncol=4, scales='free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  geom_hline(yintercept=1e-3)
g
ggsave(paste0(output.dir,'/maxGradient.png'), plot=g, height=5, width=8, units='in', dpi=250)

g <- ggplot(aic.df, aes(x=RhoConfig, y=Converge, color=Knots)) +
  theme_gray() +
  geom_point() +
  # facet_grid(Est_Spatial_RE~Species, scales='free') +
  facet_wrap(~Species, ncol=4) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
g
ggsave(paste0(output.dir,'/Convergence.png'), plot=g, height=5, width=8, units='in', dpi=250)




#Plotting Rockfish ===========================================================

heights <- 9
widths <- 7

rockfish <- c('Pacific ocean perch','Northern rockfish','Harlequin rockfish')

t <- 1
for(t in 1:n.trial.knots) {
  temp.knots <- trial.knots[t]
  
  #Estimate Spatial RE
  do.spatial <- TRUE
  temp.df <- output.df[which(output.df$Species %in% rockfish),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g
  ggsave(paste0(output.dir,'/Rockfish ', temp.knots,'kt-STRE.png'), plot=g, 
         height=heights/2, width=widths, units='in', dpi=1e3)
  
  #Don't Estimate Spatial RE
  do.spatial <- FALSE
  temp.df <- output.df[which(output.df$Species %in% rockfish),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g2 <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g2
  
  #Bring the figures together
  g12 <- plot_grid(g, g2, nrow=2, ncol=1, align='v')
  # g12
  ggsave(paste0(output.dir,'/Rockfish ', temp.knots,'kt.png'), plot=g12, 
         height=heights, width=widths, units='in', dpi=1e3)
}#next t


#Plotting Pollock and Pcod ===========================================================
t <- 1
for(t in 1:n.trial.knots) {
  temp.knots <- trial.knots[t]
  
  #Estimate Spatial RE
  do.spatial <- TRUE
  temp.df <- output.df[which(output.df$Species %in% c('Walleye pollock','Pacific cod')),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g
  ggsave(paste0(output.dir,'/Pollock_Cod ', temp.knots,'kt-STRE.png'), plot=g, 
         height=heights/2, width=widths, units='in', dpi=1e3)
  
  #Don't Estimate Spatial RE
  do.spatial <- FALSE
  temp.df <- output.df[which(output.df$Species %in% c('Walleye pollock','Pacific cod')),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g2 <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g2
  
  #Bring the figures together
  g12 <- plot_grid(g, g2, nrow=2, ncol=1, align='v')
  # g12
  ggsave(paste0(output.dir,'/Pollock_Cod ', temp.knots,'kt.png'), plot=g12, 
         height=heights, width=widths, units='in', dpi=1e3)
}#next t


#Plotting Other Species ===========================================================
t <- 1
for(t in 1:n.trial.knots) {
  temp.knots <- trial.knots[t]
  
  #Estimate Spatial RE
  do.spatial <- TRUE
  temp.df <- output.df[-which(output.df$Species %in% c('Walleye pollock','Pacific cod',rockfish)),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g
  ggsave(paste0(output.dir,'/Others ', temp.knots,'kt-STRE.png'), plot=g, 
         height=heights/2, width=widths, units='in', dpi=1e3)
  
  #Don't Estimate Spatial RE
  do.spatial <- FALSE
  temp.df <- output.df[-which(output.df$Species %in% c('Walleye pollock','Pacific cod',rockfish)),]
  temp.df <- temp.df[temp.df$Est_Spatial_RE==do.spatial | is.na(temp.df$Est_Spatial_RE),]
  temp.df$Species <- gsub(" ", "\n", temp.df$Species)
  
  g2 <- ggplot(temp.df, aes(x=Year, y=value, fill=Region)) +
    theme_dark()+
    theme(legend.position='top') +
    scale_fill_colorblind() +
    geom_area(position='stack', alpha=alpha) +
    facet_grid(Species~RhoConfig) +
    ggtitle(paste('Estimate Spatial RE:', do.spatial)) +#, 
    # subtitle=paste('Knots:',temp.knots)) +
    ylab('Proportion') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  # g2
  
  #Bring the figures together
  g12 <- plot_grid(g, g2, nrow=2, ncol=1, align='v')
  # g12
  ggsave(paste0(output.dir,'/Others ', temp.knots,'kt.png'), plot=g12, 
         height=heights, width=widths, units='in', dpi=1e3)
}#next t

#PLOT: Spatial Distribution Examples ==============================================
n.species
species

loc.rockfish <- which(species.list$name %in% rockfish)

#Extract the model fit

names(vast_est.output[[1]][[3]])

#Harlequin Rockfish
Opt <- vast_est.output[[1]][[loc.rockfish[3]]]


MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"="Gulf_of_Alaska",
                                                   "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
                                                   "Extrapolation_List"=Extrapolation_List )

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




#===========================================================
#Printing GOA POP Table
# write.csv(vast_est.output[[7]][[which(species.list$name %in% 'Pacific ocean perch')]]$vast_est, file=paste0(output.dir,'/POP AR+RW spat_TRUE.csv'))

#===========================================================
#Plotting All Species separately

# spec.surv <- unique(output.df$Name)
# n.spec.surv <- length(spec.surv)
# 
# s <- 1
# for(s in 1:n.spec.surv) {
#   temp.name <- spec.surv[s]
#   
#   temp.df <- output.df[output.df$Name==temp.name,]
#   
#   t <- 1
#   for(t in 1:n.trial.knots) {
#     temp.knots <- trial.knots[t]
#     g <- ggplot(temp.df[temp.df$Knots==temp.knots | temp.df$Knots==FALSE,],
#                   aes(x=Year, y=value, fill=Region)) +
#            theme_dark() +
#            theme(legend.position='top') +
#            scale_fill_colorblind() +
#            ylab('Proportion') +
#            geom_area(position='stack', alpha=alpha) +
#            facet_grid(Est_Spatial_RE~RhoConfig) +
#            ggtitle(temp.name, subtitle=paste('Knots:',temp.knots)) 
#     g
#   }#next t
# }#next s

