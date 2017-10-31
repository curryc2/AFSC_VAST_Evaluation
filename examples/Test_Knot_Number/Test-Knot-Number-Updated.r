#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Compare Knot Number --- UPDATED FOR LOWER STORAGE REQUIREMENTS
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 4.10.17
#
#Purpose: To explore sensitivity of model-based index estimates to different knot number specification
#             THIS SCRIPT IS AN UPDATE FROM EARLIER VERSION TO REDUCE STORAGE REQUIREMENTS.
#
#
#
#==================================================================================================
#NOTES:
#
#  a) Could calculate design-based estiamtes within wrapper function or at the end while plotting.
#  b) Trouble with bias.corr=TRUE, Memory allocation fail, 200 kts, unknown species, during parallel. 
#  b.1) Problem is with EBS Arrowtooth. for bias.corr=TRUE, this species will be remove and knot numbers limited. 
#==================================================================================================
#TIMING:
# For 10-1000, by 100 knots
# [1] "### START: Mon Apr 10 16:17:08 2017"
# [1] "### END: Mon Apr 10 22:54:04 2017"
#==================================================================================================

require(VAST)
require(TMB)
require(parallel)
require(snowfall)
require(ggplot2)
require(cowplot)
require(xlsx)
require(ggthemes)


source("R/calc-design-based-index.r")
source("R/create-VAST-input.r")
source("R/create-Data-Geostat.r")
source("R/load-RACE-data.r")

source("R/cleanup-VAST-file.r")
source("R/get-VAST-index.r")


home.dir <- getwd()
#Create working directory
working.dir <- paste0(home.dir, "/examples/Test_Knot_Number")

#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)

#Limit species included
species.list <- species.list[species.list$include=='Y',]
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
trial.knots <- c(100,200,300,400,500,750,1000)#seq(100, 1000, by=100)
n.trial.knots <- length(trial.knots)

#Boolean for bias correction
bias.correct <- FALSE

#Update if
if(bias.correct==TRUE) {
  species.list <- species.list[-which(species.list$survey=='EBS_SHELF'),]
  n.species <- nrow(species.list)
  species.series <- c(1:n.species)
  
  
  trial.knots <- c(100,200,300,400,500)
  n.trial.knots <- length(trial.knots)
}
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
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Delta1 = 0, Delta2 = 0)

ObsModel = c(1, 0) #Lognormal

#SPECIFY OUTPUTS
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
            Calculate_Cov_SE = 0, Calculate_Synchrony = 0,
            Calculate_Coherence = 0)

#Output Directory Name
output.dir <- paste0(working.dir,"/output_bias.correct_",bias.correct)



#=======================================================================
##### WRAPPER FUNCTION FOR RUNNING IN PARALLEL #####


s <- 1
# for(s in 1:n.species) {
species_wrapper_fxn_knots <- function(s, n_x) {
  
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
     "TmbList", "Obj", "Opt", "Report", "Save")
  
  #========================================================================
  setwd(home.dir)
  ##### RETURN SECTION #####
  return(vast_est)
} 


#=======================================================================
##### Loop Through Trial Knots  #####
if(do.estim==TRUE) {
  vast_est.output <- vector('list', length=n.trial.knots)

  time.1 <- date()
  
  t <- 1
  for(t in 1:n.trial.knots) {
    print(paste('## Trial Knot Number',t,'of',n.trial.knots))
    print(paste('# Trial Knots:',trial.knots[t]))
    #Specify trial observation model
    
    #Specify knots
    n_x <- trial.knots[t]
    
    #Setup File
    trial.dir <- paste0(working.dir,"/",n_x,"_bias.corr_",bias.correct)
    dir.create(trial.dir)
    
    
    #=======================================================================
    ##### SNOWFALL CODE FOR PARALLEL #####
    sfInit(parallel=TRUE, cpus=n.cores, type='SOCK')
    sfExportAll() #Exportas all global variables to cores
    sfLibrary(TMB)  #Loads a package on all nodes
    sfLibrary(VAST)
    output <- sfLapply(species.series, fun=species_wrapper_fxn_knots, n_x=n_x)
    # sfRemove(Save)
    # sfRemover(VAST_input)
    sfStop()
    
    vast_est.output[[t]] <- output
    
  }# next t
  
  #Dimensions for vast_est.output are 1) Trial knots, 2) Species
  # vast_est.output[[1:n.trial.knots]][[1:n.species]]
  
  #Create output directory
  dir.create(output.dir)
  save(vast_est.output, file=paste0(output.dir,"/vast_est.output.RData"))
  
  #=======================================================================
  ##### DELETE UNNECESSARY FILE STRUCTURE #####
  #Must reset working directory
  setwd(working.dir)
  t <- 1
  for(t in 1:n.trial.knots) {
    # s <- 1
    # for(s in 1:n.species) {
      # temp.DateFile <- paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct,"/",species.list$name[s],"/")
      #Remove manually
      # file.remove(paste0(temp.DateFile,"parameter_estimates.Rdata"))
      # file.remove(paste0(temp.DateFile,"parameter_estimates.txt"))
      
      # #REMOVE EVERYTHING...
      # # file.remove(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct,"/"), force=TRUE)
      # print(unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct,"/"), recursive=TRUE))
      
      
      # unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct,"/",species.list$name[s]), recursive=TRUE)
    # }#next s
    
    unlink(paste0(working.dir,"/",trial.knots[t],"_bias.corr_",bias.correct), recursive=TRUE)
  }#next t
  
  time.2 <- date()
  
  print(paste('### START:', time.1))
  print(paste('### END:', time.2))
  
}else {
  load(paste0(output.dir,"/vast_est.output.RData"))
}

#=======================================================================
##### Create Output Lists #####

#Create large list
vast.list  <- NULL
#Year, Biomass, SD, CV, Species, Model, Knots

t <- 1
for(t in 1:n.trial.knots) {
  s <- 1
  for(s in 1:n.species) {
    
    temp.list <- vast_est.output[[t]][[s]][c(1,4,6)]
    #Calculate CV
    CV <- temp.list$SD_mt/temp.list$Estimate_metric_tons
    temp.species <- species.list$name[s]
    temp.survey <- species.list$survey[s]
    temp.name <- paste0(temp.survey,": ",temp.species)
    #Bind it
    temp.list <- cbind(temp.list, CV, temp.survey, temp.species, temp.name, 'VAST', trial.knots[t])
    #Add it
    vast.list <- rbind(vast.list, temp.list)
  }#next s
}#next t
names(vast.list) <- c('Year','Biomass','SD','CV','Survey','Species','Name','Model','Knots')


#Add Design-based estimates
db.list <- NULL
s <- 1
for(s in 1:n.species ) {
  #Get design-based estimate
  db_est <- calc_design_based_index(species.codes=species.list$species.code[s], survey=species.list$survey[s])
  temp.species <- species.list$name[s]
  temp.survey <- species.list$survey[s]
  temp.name <- paste0(temp.survey,": ",temp.species)
  temp.list <- cbind(db_est[,c(1,2,4,5)], temp.survey, temp.species, temp.name, "Design-based", 'Design-based')
  #Add it
  db.list <- rbind(db.list, temp.list)
}#next s
names(db.list) <- c('Year','Biomass','SD','CV','Survey','Species', 'Name','Model','Knots')

#Combine the lists
survey.list <- rbind(vast.list,db.list)


#=======================================================================
##### Plot Comparison of Results #####
scale.hues <- c(5,250)
dpi <- 500
#Find years survey was conducted

surveys <- unique(survey.list$Survey)
n.surveys <- length(surveys)

survey.list$Knots <- factor(survey.list$Knots, levels=c(trial.knots,'Design-based'))

#========================================
#Plot All Species Across Surveys
s <- 1
for(s in 1:n.surveys) {
#GOA
survey <- surveys[s]
#Limit data set
plot.list <- survey.list[survey.list$Survey==survey,]
yrs.surv <- sort(unique(plot.list$Year[plot.list$Model=='Design-based']))
plot.list <- plot.list[plot.list$Year %in% yrs.surv,]

#Remove 2001 from design-based results because of incomplete sampling
if(survey=='GOA') { plot.list <- plot.list[-which(plot.list$Year==2001 & plot.list$Model=='Design-based'),] }
plot.list$Biomass <- plot.list$Biomass/1e3

#PLOT Indices
g <- ggplot(plot.list, aes(x=Year, y=Biomass, color=Knots, lty=Model, ymin=0)) +
       theme_gray() +
       # theme_economist() +
       theme(legend.position='right') +
       geom_line() +
       facet_wrap(~Species, scales='free', ncol=2) +
       labs(list(y='Biomass (thousands of metric tonnes)')) +
       # ggtitle('Survey:', subtitle='Gulf of Alaska') +
       ggtitle(paste(survey, 'Survey')) +
       scale_color_hue(h=scale.hues)
       # scale_color_brewer(type='seq', palette=1)
g

 

# g
ggsave(paste0(output.dir,"/", survey," VAST Index Compare v DB.png"), g, height=9, width=8, units='in', dpi=dpi)

#Vast Models only
g2 <- ggplot(plot.list[plot.list$Model=='VAST',], aes(x=Year, y=Biomass, color=Knots, ymin=0)) +
        theme_gray() +
        geom_line() +
        theme(legend.position='right') +
        facet_wrap(~Species, scales='free', ncol=2) +
        labs(list(y='Biomass (thousands of metric tonnes)')) +
        # ggtitle('Survey:', subtitle='Gulf of Alaska') +
        ggtitle(paste(survey, 'Survey')) +
        scale_color_hue(h=scale.hues)

# g2
ggsave(paste0(output.dir,"/", survey," VAST Index Compare.png"), g2, height=9, width=8, units='in', dpi=dpi)

#Plot Survey Variance Measures

g3 <- ggplot(plot.list, aes(x=Species, y=CV, fill=Knots)) +
        theme_gray() +
        theme(legend.position='right') +
        # geom_boxplot(aes(lty=Model)) +
        geom_boxplot(aes(color=Model)) +
        scale_color_colorblind() +
        # scale_color_hue() +
        labs(list(y=paste('Annual Survey CV'))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, debug=FALSE)) +
        scale_fill_hue(h=scale.hues) +
        ggtitle(paste(survey, 'Survey'))
  
        
g3
ggsave(paste0(output.dir,"/", survey," CV Compare.png"), g3, height=5, width=7, units='in', dpi=dpi)
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/", survey," CV Compare.png"), g3, height=5, width=8, units='in', dpi=1e3)

#Separate boxes
g4 <- ggplot(plot.list, aes(x=Species, y=CV, fill=Knots)) +
  theme_gray() +
  theme(legend.position='bottom') +
  # geom_boxplot(aes(lty=Model)) +
  geom_boxplot(aes(color=Model)) +
  scale_color_colorblind() +
  # scale_color_hue() +
  labs(list(y=paste('Annual Survey CV'))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_hue(h=scale.hues) +
  # ggtitle(paste(survey, 'Survey')) +
  facet_wrap(~Species, scales='free', nrow=2)
  


g4
ggsave(paste0(output.dir,"/", survey," CV Compare.png"), g4, height=6, width=8, units='in', dpi=dpi)
ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/", survey, " CV Compare_2.png"), g4, height=7, width=8, units='in', dpi=1000)

#Plotting for GPT presentation
g5 <- ggplot(plot.list, aes(x=Year, y=Biomass, color=Knots, lty=Model, ymin=0)) +
        theme_gray() +
        # theme_wsj() +
        # theme_solarized() +
        theme(legend.position='right') +
        geom_line() +
        facet_wrap(~Species, scales='free', ncol=3) +
        labs(list(y='Biomass (thousands of metric tonnes)')) +
        # ggtitle('Survey:', subtitle='Gulf of Alaska') +
        ggtitle(paste(survey, 'Survey')) +
        scale_color_hue(h=scale.hues) +
        geom_line(data=plot.list[plot.list$Model=='Design-based',], color='black') +
        geom_point(data=plot.list[plot.list$Model=='Design-based',], show.legend=FALSE, colour='black')

g5

ggsave(paste0(getwd(),"/Output/Figs for Sept_2017 GPT/", survey, " Index Compare.png"), g5, height=7, width=10, units='in', dpi=1000)
# png(paste0(output.dir,'/', survey, ' Figures.png'), height=7, width=8, units='in', res=500)
# dev.off()
}

#==================================================================================
#GENERATE TABLES OF % DIFFERENCE IN CV'S BETWEEN METHODS

#Gulf of Alaska
temp.survey <- 'GOA'
temp.db <- db.list[db.list$Survey==temp.survey,]
temp.vast <- vast.list[vast.list$Survey==temp.survey,]

temp.species <- unique(temp.db$Species)
n.temp.species <- length(temp.species)

temp.years <- sort(unique(temp.db$Year))
n.temp.years <- length(temp.years)

pct.diff.cv.goa <- array(dim=c(n.temp.species,n.trial.knots,n.temp.years), dimnames=list(temp.species,trial.knots,temp.years))
pct.diff.bio.goa <- array(dim=c(n.temp.species,n.trial.knots,n.temp.years), dimnames=list(temp.species,trial.knots,temp.years))
s <- 1
for(s in 1:n.temp.species) {
  t <- 1
  for(t in 1:n.trial.knots) {
    y <- 1
    for(y in 1:n.temp.years) {
      #CV
      pct.diff.cv.goa[s,t,y] <- (temp.vast$CV[temp.vast$Species==temp.species[s] & temp.vast$Knots==trial.knots[t] & temp.vast$Year==temp.years[y]] - 
                                   temp.db$CV[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]])/
                                   temp.db$CV[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]]
      #Biomass Estimate
      pct.diff.bio.goa[s,t,y] <- (temp.vast$Biomass[temp.vast$Species==temp.species[s] & temp.vast$Knots==trial.knots[t] & temp.vast$Year==temp.years[y]] - 
                                    temp.db$Biomass[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]])/
                                    temp.db$Biomass[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]]
    }#next y
  }#next t
}#next s

#Gulf of Alaska
temp.survey <- 'AI'
temp.db <- db.list[db.list$Survey==temp.survey,]
temp.vast <- vast.list[vast.list$Survey==temp.survey,]

temp.species <- unique(temp.db$Species)
n.temp.species <- length(temp.species)

temp.years <- sort(unique(temp.db$Year))
n.temp.years <- length(temp.years)

pct.diff.cv.ai <- array(dim=c(n.temp.species,n.trial.knots,n.temp.years), dimnames=list(temp.species,trial.knots,temp.years))
pct.diff.bio.ai <- array(dim=c(n.temp.species,n.trial.knots,n.temp.years), dimnames=list(temp.species,trial.knots,temp.years))
s <- 1
for(s in 1:n.temp.species) {
  t <- 1
  for(t in 1:n.trial.knots) {
    y <- 1
    for(y in 1:n.temp.years) {
      #CV
      pct.diff.cv.ai[s,t,y] <- (temp.vast$CV[temp.vast$Species==temp.species[s] & temp.vast$Knots==trial.knots[t] & temp.vast$Year==temp.years[y]] - 
                                  temp.db$CV[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]])/
                                  temp.db$CV[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]]
      #Biomass Estimate
      pct.diff.bio.ai[s,t,y] <- (temp.vast$Biomass[temp.vast$Species==temp.species[s] & temp.vast$Knots==trial.knots[t] & temp.vast$Year==temp.years[y]] - 
                                   temp.db$Biomass[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]])/
                                   temp.db$Biomass[temp.db$Species==temp.species[s] & temp.db$Year==temp.years[y]]
    }#next y
  }#next t
}#next s

#PRINT TABLE

write.csv(apply(pct.diff.cv.goa, c(1,2), mean), file=paste0(output.dir, "/GOA Mean Pct Diff CV.csv"))
write.csv(apply(pct.diff.cv.ai, c(1,2), mean), file=paste0(output.dir, "/AI Mean Pct Diff CV.csv"))

#==================================================================================
#Plot GOA Pollock
survey <- 'Gulf of Alaska'
plot.list <- survey.list[survey.list$Survey=='GOA' & survey.list$Species=='Walleye pollock',]
yrs.surv <- sort(unique(plot.list$Year[plot.list$Model=='Design-based']))
plot.list <- plot.list[plot.list$Year %in% yrs.surv,]

#Remove 2001 from design-based results because of incomplete sampling
if(survey=='GOA') { plot.list <- plot.list[-which(plot.list$Year==2001 & plot.list$Model=='Design-based'),] }
plot.list$Biomass <- plot.list$Biomass/1e3

#PLOT Indices
g.idx <- ggplot(plot.list, aes(x=Year, y=Biomass, color=Knots, lty=Model)) +
           theme_gray() +
            # theme(legend.position='bottom') +
           geom_line() +
           facet_wrap(~Species, scales='free') +
           labs(list(y='Biomass (thousands of metric tonnes)')) +
           # ggtitle('Survey:', subtitle='Gulf of Alaska') +
           # ggtitle(paste(survey, 'Survey')) +
           scale_color_hue(h=scale.hues)
g.idx
ggsave(paste0(output.dir,"/GOA Pollock Idx.png"), g.idx, height=6, width=8, units='in', dpi=dpi)

g.cv <- ggplot(plot.list, aes(x=Knots, y=CV, fill=Knots)) +
  theme_gray() +
  # geom_boxplot(aes(lty=Model)) +
  geom_boxplot() +
  facet_wrap(~Species, scales='free') +
  labs(list(y=paste('Annual Survey CV'))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, debug=FALSE)) +
  scale_fill_hue(h=scale.hues) +
  # ggtitle(paste(survey, 'Survey')) +
  
  theme(legend.position='none')
  
# g.cv


#Combine plots with gridExtra
g.both <- plot_grid(g.idx, g.cv, nrow=1, ncol=2, rel_widths=c(3.5,1))
g.both

ggsave(paste0(output.dir,"/GOA Pollock Idx and CV.png"), g.both, height=5, width=8, units='in', dpi=dpi)


#==================================================================================
#Plot GOA POP
survey <- 'Gulf of Alaska'
plot.list <- survey.list[survey.list$Survey=='GOA' & survey.list$Species=='Pacific ocean perch',]
yrs.surv <- sort(unique(plot.list$Year[plot.list$Model=='Design-based']))
plot.list <- plot.list[plot.list$Year %in% yrs.surv,]

#Remove 2001 from design-based results because of incomplete sampling
if(survey=='GOA') { plot.list <- plot.list[-which(plot.list$Year==2001 & plot.list$Model=='Design-based'),] }
plot.list$Biomass <- plot.list$Biomass/1e3

#PLOT Indices
g.idx <- ggplot(plot.list, aes(x=Year, y=Biomass, color=Knots, lty=Model)) +
  theme_gray() +
  # theme(legend.position='bottom') +
  geom_line() +
  facet_wrap(~Species, scales='free') +
  labs(list(y='Biomass (thousands of metric tonnes)')) +
  # ggtitle('Survey:', subtitle='Gulf of Alaska') +
  # ggtitle(paste(survey, 'Survey')) +
  scale_color_hue(h=scale.hues)
g.idx
ggsave(paste0(output.dir,"/GOA POP Idx.png"), g.idx, height=6, width=8, units='in', dpi=dpi)

g.cv <- ggplot(plot.list, aes(x=Knots, y=CV, fill=Knots)) +
  theme_gray() +
  # geom_boxplot(aes(lty=Model)) +
  geom_boxplot() +
  facet_wrap(~Species, scales='free') +
  labs(list(y=paste('Annual Survey CV'))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, debug=FALSE)) +
  scale_fill_hue(h=scale.hues) +
  # ggtitle(paste(survey, 'Survey')) +
  
  theme(legend.position='none')

# g.cv


#Combine plots with gridExtra
g.both <- plot_grid(g.idx, g.cv, nrow=1, ncol=2, rel_widths=c(3.5,1))
g.both

ggsave(paste0(output.dir,"/GOA POP Idx and CV.png"), g.both, height=5, width=8, units='in', dpi=dpi)

#GENERATE WORKSHEET FOR WITH ALL OF THE GOA POP INDICES
write.xlsx(survey.list[survey.list$Survey=='GOA' & survey.list$Species=='Pacific ocean perch' &
                         survey.list$Year %in% yrs.surv,],
             file=paste0(output.dir,'/GOA POP Indices.xlsx'), sheetName='plot.list')



#==================================================================================
#Plot GOA Arrowtooth for Ingrid's Assessment
survey <- 'Gulf of Alaska'
plot.list <- survey.list[survey.list$Survey=='GOA' & survey.list$Species=='Arrowtooth flounder',]
yrs.surv <- sort(unique(plot.list$Year[plot.list$Model=='Design-based']))
plot.list <- plot.list[plot.list$Year %in% yrs.surv,]

#Remove 2001 from design-based results because of incomplete sampling
plot.list <- plot.list[-which(plot.list$Year==2001 & plot.list$Model=='Design-based'),]
plot.list <- plot.list[-which(plot.list$Year==2017 & plot.list$Model=='Design-based'),]

plot.list$Biomass <- plot.list$Biomass/1e3

#PLOT Indices
g.idx <- ggplot(plot.list, aes(x=Year, y=Biomass, color=Knots, lty=Model)) +
  # theme_gray() +
  # theme(legend.position='bottom') +
  # geom_line() +
  # facet_wrap(~Species, scales='free') +
  # labs(list(y='Biomass (thousands of metric tonnes)')) +
  # # ggtitle('Survey:', subtitle='Gulf of Alaska') +
  # # ggtitle(paste(survey, 'Survey')) +
  # scale_color_hue(h=scale.hues)




  theme_gray() +
  # theme_wsj() +
  # theme_solarized() +
  geom_line() +
  facet_wrap(~Species, scales='free', ncol=3) +
  labs(list(y='Biomass (thousands of metric tonnes)')) +
  # ggtitle('Survey:', subtitle='Gulf of Alaska') +
  ggtitle(paste(survey, 'Survey')) +
  scale_color_hue(h=scale.hues) +
  geom_line(data=plot.list[plot.list$Model=='Design-based',], color='black') +
  geom_point(data=plot.list[plot.list$Model=='Design-based',], show.legend=FALSE, colour='black')

# g.idx

g.cv <- ggplot(plot.list, aes(x=Knots, y=CV, fill=Knots)) +
  theme_gray() +
  # geom_boxplot(aes(lty=Model)) +
  geom_boxplot() +
  facet_wrap(~Species, scales='free') +
  labs(list(y=paste('Annual Survey CV'))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, debug=FALSE)) +
  scale_fill_hue(h=scale.hues) +
  # ggtitle(paste(survey, 'Survey')) +
  
  theme(legend.position='none')

  # theme_gray() +
  # theme(legend.position='right') +
  # # geom_boxplot(aes(lty=Model)) +
  # geom_boxplot(aes(color=Model)) +
  # scale_color_colorblind() +
  # # scale_color_hue() +
  # labs(list(y=paste('Annual Survey CV'))) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, debug=FALSE)) +
  # scale_fill_hue(h=scale.hues) +
  # ggtitle(paste(survey, 'Survey'))
# g.cv


#Combine plots with gridExtra
g.both <- plot_grid(g.idx, g.cv, nrow=1, ncol=2, rel_widths=c(3.5,1))
g.both

ggsave(paste0(output.dir,"/GOA Arrowtooth Index and CV.png"), g.both, height=5, width=8, units='in', dpi=dpi)


#Facet
# g.multi <- vector('list', length=n.species)
# 
# s <- 1
# for(s in 1:n.species) {
#   g.multi[[s]] <- ggplot(plot.list[plot.list$Species=='Big skate',], aes(x=Species, y=CV, fill=Knots)) +
#     theme_gray() +
#     geom_boxplot() +
#     labs(list(y='Annual Survey CV')) +
#     facet_wrap(~Species, ncol=5, drop=TRUE)
#     
#     if(s > 1) { +  }
#     + theme(legend.position="none")
# }
# 
# 
# require(gridExtra)
# 
# do.call("grid.arrange", c(g.multi, ncol=5))

# 
# ggsave(paste0(output.dir, "/VAST Index Compare.png"), g2, height=6, width=9, units='in')
# 







