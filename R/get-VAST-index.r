#' @title
#' Get VAST model-based survey index and uncertainty
#' 
#' @description 
#' \code{get_VAST_index} returns table of index values and errors.
#'
#' @param TmbData Input data object for VAST model
#' @param Sdreport Output object with parameter and uncertainty (delta-method) estimates from vAST
#' @param bias.correct Boolean for whether bias corrected results should be returned
#' @param Data_Geostat 
#'
#' @return A table of index values, by year, with associated standard deviation

#' @export
#'
get_VAST_index <- function(TmbData, Sdreport=Opt[["SD"]], bias.correct, Data_Geostat=NULL) {
  ###Testing
  Sdreport=Opt[["SD"]]
  bias.correct <- FALSE
  
  #Initial Checkup
  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    ParName = "Index_tl"
    TmbData[['n_c']] = 1
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    ParName = "Index_ctl"
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    ParName = "Index_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }
  if( "Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialVAM
    ParName = "Index_tp"
    TmbData[["n_l"]] = 1
    TmbData[["n_c"]] = TmbData[["n_p"]]
  }
  
  #Retreive Years
  if(is.null(Data_Geostat)) {
    Year_Set = 1:TmbData$n_t
    Years2Include = 1:TmbData$n_t
  }else {
    Year_Set <- seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
    Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  }
  strata_names = 1:TmbData$n_l
  category_names = 1:TmbData$n_c
  
  # Extract index
  if( ParName %in% c("Index_tl","Index_ctl","Index_cyl")){
    if( bias.correct==TRUE && "unbiased"%in%names(Sdreport) ){
      log_Index_ctl = array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==paste0("ln_",ParName))],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),'Std. Error']), dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
      Index_ctl = array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==ParName)],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),'Std. Error']), dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
    }else{
      log_Index_ctl = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),], dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
      Index_ctl = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),], dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
    }
  }
  if( ParName %in% c("Index_tp")){
    if( bias.correct==TRUE && "unbiased"%in%names(Sdreport) ){
      Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==ParName)],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==paste0("ln_",ParName))],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'SpatialVAM'" )
      }
    }else{
      Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'SpatialVAM'" )
      }
    }
  }
  
  Table = NULL
  for( cI in 1:TmbData$n_c ){
    Tmp = data.frame( "Year"=Year_Set, "Unit"=1, "Fleet"=rep(strata_names,each=TmbData$n_t), "Estimate_metric_tons"=as.vector(Index_ctl[cI,,,'Estimate']), "SD_log"=as.vector(log_Index_ctl[cI,,,'Std. Error']), "SD_mt"=as.vector(Index_ctl[cI,,,'Std. Error']) )
    if( TmbData$n_c>1 ) Tmp = cbind( "Category"=category_names[cI], Tmp)
    Table = rbind( Table, Tmp )
  }
  
  
  return(Table)
  
  
}