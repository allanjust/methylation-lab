#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  

#' recall that we have a processed dataset with 15 samples
dim(WB)

suppressPackageStartupMessages({
  library(CpGassoc)
})

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on sex (as a predictor of methylation)

#' and with adjustment for cell types

#' looking at results


#' End of script