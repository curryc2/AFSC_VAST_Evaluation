run_RE_model <- function(input.yrs, input.idx, input.cv, DateFile, home.dir, n_PE=1, PE_vec=c(1,1,1)) {
  require(R2admb)
  require(PBSmodelling)
  ### TESTING ###
  # n_PE <- 3#1
  # PE_vec <- 1:3#c(1,1,1)
  ###############
  
  setwd(DateFile)
  
  ### CREATE ADMB INPUT FILE ###
  styr <- min(input.yrs)  #Start year
  endyr <- max(input.yrs)  #End year
  num_indx <- ncol(input.idx) #Number of survey indices
  # n_PE <- 1 #3  #Number of process error parameters
  nobs <- length(input.yrs)  #Number of surveys
  yrs_srv <- input.yrs  #Survey years
  srv_est <- input.idx  #Survey Biomass
  srv_cv <- input.cv  #Survey Biomass CV

  write_dat(paste0(DateFile,'/re'), L=list(styr, endyr, num_indx, n_PE, PE_vec,
                               nobs, yrs_srv,
                               srv_est, srv_cv))
  
  ### COPY IN .tpl ###
  file.copy(from=paste0(home.dir,"/admb/re.tpl"), to=paste0(DateFile,"re.tpl"), overwrite=TRUE)
  
  ### COMPILE ADMB-RE MODEL ###
  compile_admb(fn="re", re=TRUE, verbose=FALSE)
  
  ### ADMB RANDOM EFFECTS MODEL ###
  run_admb(fn="re", verbose=FALSE)
  
  ### READ MODEL OUTPUTS ###
  biomA <- readList(paste0(DateFile,"/biomA.out"))$biomA
  
  return(biomA)
}