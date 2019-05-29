#'## Install packages needed for epigenomics tutorial  
#' This should be run before the workshop
#' because of the need to download large datasets.  
#'   

#'## It is important that you have already updated R
#'you should be running version 3.6.0
R.version$version.string   
#' this is because many packages update and change to fix bugs and add new features.
#' Local library   
.libPaths("C:/EBC4/Rpackages")

#' Please consider updating your packages
#' this step requires agreeing for each update (y for yes)
update.packages(lib.loc="C:/EBC4/Rpackages")

#'# Installation of new packages   
#' vector of packages we will need if not yet installed:

methpackagesCRAN <- c("Rcpp","openssl","CpGassoc", "rmarkdown", "knitr", "matrixStats","reshape","glmnet","statmod","XML","BiocManager",
                      "pryr", "data.table", "qqman", "RPMM", "MASS", "sandwich", "lmtest","foreach", "stringi","doParallel","magrittr","purrr")
methpackagesBioC <- c("minfi", "FlowSorted.Blood.450k", "missMethyl", "ENmix","IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICmanifest",
                      "sva", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", 
                      "DMRcate", "shinyMethyl","bumphunter","wateRmelon","FDb.InfiniumMethylation.hg19")
#' install these from CRAN:
toinstallCRAN <- setdiff(methpackagesCRAN, installed.packages()[,1])
if(length(toinstallCRAN >= 1)) {
  install.packages(toinstallCRAN,dependencies=TRUE, lib="C:/EBC4/Rpackages")
  cat("finished installing new packages from CRAN\n")
} else cat("packages we need from CRAN are already installed\n")
#' install these from BioConductor:
toinstallBioC <- setdiff(methpackagesBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  BiocManager::install(toinstallBioC, update = TRUE ,lib.loc="C:/EBC4/Rpackages")
  cat("finished installing new packages from BioConductor\n")
} else cat("packages we need from BioConductor are already installed\n")

devtools::install_github("hhhh5/ewastools@v1.5")

#' check that we were successful
if(!all(c(toinstallBioC, toinstallCRAN) %in% installed.packages()[,1])) stop(
  "required packages not installed - please retry script carefully making sure you have already updated R and work through any error messages")

if(!as.numeric(sub("\\.[0-9]$", "", installed.packages()["minfi","Version"])) >= 1.24) stop(
  "you don't have the minfi version needed for this workshop")

#' Session Information
#' If you cannot successfully work through this script, please run the following two commands 
#' and send the output to the workshop organizers with your request for help:
#sessionInfo()
#installed.packages()[,c("Package", "Version")]

#' cleanup
rm(methpackagesCRAN, methpackagesBioC, toinstallCRAN, toinstallBioC)
#' End of script

