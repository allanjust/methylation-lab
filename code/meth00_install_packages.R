#'## Install packages needed for epigenomics tutorial  
#' This should be run before the workshop
#' because of the need to download large datasets.  
#'   

#'## It is important that you have already updated R
#'you should be running version 3.3.1
R.version$version.string   
#' this is because many packages update and change to fix bugs and add new features.  
#'   
#' vector of packages we will need if not yet installed:
methpackagesCRAN <- c("CpGassoc", "ggplot2", "matrixStats", "pryr")
methpackagesBioC <- c("minfi", "FlowSorted.CordBlood.450k", "missMethyl", "LOLA", "coMET")
#' install these from CRAN:
toinstallCRAN <- setdiff(methpackagesCRAN, installed.packages()[,1])
if(length(toinstallCRAN >= 1)) install.packages(toinstallCRAN)
#' install these from BioConductor:
toinstallBioC <- setdiff(methpackagesBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(toinstallBioC)
}  
#' cleanup
rm(methpackagesCRAN, methpackagesBioC, toinstallCRAN, toinstallBioC)

#' Session Information
sessionInfo()

#' End of script