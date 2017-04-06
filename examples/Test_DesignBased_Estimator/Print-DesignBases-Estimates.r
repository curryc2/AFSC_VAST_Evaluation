#==================================================================================================
#Project Name: VAST spatial delta-GLMM (Thorson) Evaluation: Print Design-Based Estimates to Ensure Correct
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 4.6.17
#
#Purpose: To save design-based estimates to .csv for sharing with assessment authors for confirmation.. OK
#
#
#
#==================================================================================================
#NOTES:
#
#
#==================================================================================================
require(xlsx)
source("R/calc-design-based-index.r")

Region <- "Gulf_of_Alaska"

#Determine species list
species.list <- read.csv("data/eval_species_list.csv", stringsAsFactors=FALSE)
#Limit to those included
species.list <- species.list[species.list$include=='Y',]

n.species <- nrow(species.list)
species.series <- c(1:n.species)



s <- 1
for(s in 1:n.species) {
  print(paste(s, 'of', n.species))
  species.codes <- species.list$species.code[s]
  
  spec.name <- species.list$name[s]
  db.est <- calc_design_based_index(species.codes=species.codes, Region=Region)
  write.xlsx(db.est, file=paste0(getwd(), "/examples/Test_DesignBased_Estimator/Design Based Estimates.xlsx"), sheetName=spec.name,
               append=ifelse(s==1,FALSE,TRUE))
  
}
  
