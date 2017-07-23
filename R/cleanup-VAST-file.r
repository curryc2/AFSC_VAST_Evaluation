#' Delete extra VAST files to reduce overhead, after estimation has been completed 
#'
#' @param DateFile string with file extension for current directory
#' @param Version string identifying VAST version
#'
#' @export
cleanup_VAST_file <- function(DateFile, Version="VAST_v2_4_0") {
  #Delete VAST .o ~ 14,387 KB
  unlink(paste0(DateFile,Version,".o"))
  unlink(paste0(DateFile,Version,".cpp"))
  
  
  #Delete VAST .dll ~ 13,3559 KB
  if(.Platform$OS.type=='windows') {
    unlink(paste0(DateFile,Version,".dll"))
  }
  if(.Platform$OS.type=='unix') {
    unlink(paste0(DateFile,Version,".so"))
  }
}