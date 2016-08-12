#'## Install packages needed for epigenomics tutorial  
#' This should be run before the workshop
#' because of the need to download large datasets.  
#'   

#'## It is important that you have already updated R
#'you should be running version 3.3.1
R.version$version.string   
#' this is because many packages update and change to fix bugs and add new features.  

#' Please consider updating your packages
#' this step requires agreeing for each update (y for yes)
update.packages()

#'# Installation of new packages   
#' vector of packages we will need if not yet installed:
methpackagesCRAN <- c("CpGassoc", "rmarkdown", "knitr", "ggplot2", "matrixStats", 
                      "pryr", "data.table","qqman", "RPMM")
methpackagesBioC <- c("minfi", "FlowSorted.CordBlood.450k", "missMethyl", "LOLA", "coMET", "ENmix",
                      "sva", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "DMRcate")
#' install these from CRAN:
toinstallCRAN <- setdiff(methpackagesCRAN, installed.packages()[,1])
if(length(toinstallCRAN >= 1)) {
  install.packages(toinstallCRAN)
  cat("finished installing new packages from CRAN\n")
} else cat("packages we need from CRAN are already installed\n")
#' install these from BioConductor:
toinstallBioC <- setdiff(methpackagesBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(toinstallBioC, suppressUpdates = T)
  cat("finished installing new packages from BioConductor\n")
} else cat("packages we need from BioConductor are already installed\n")
#' cleanup
rm(methpackagesCRAN, methpackagesBioC, toinstallCRAN, toinstallBioC)

#' Session Information
sessionInfo()

#' End of script