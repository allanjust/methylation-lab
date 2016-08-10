#'# Interpret methylation results  
#' Using results from CpGassoc generated in our script:  
#'  meth02_analyze_data.R  

suppressPackageStartupMessages({
  library(missMethyl) # pathway analyses
  library(LOLA) # enrichment analyses
})

#' recall that we have already run an Epigenome-Wide Association Study
#' 
nrow(results2$results[results2$results$FDR < 0.5, ])

#'# Pathway analyses with missMethyl
#' we need to select a set of CpGs that we think may have signal
#' for example, we select those with FDR < 0.1 and |delta beta| > 0.05


#' we use the gometh function from the missMethyl package
#' with a list of all of the tested CpG sites (background)
#' to account for the array design (more CpGs put in some genes)


#'# Enrichment analyses with LOLA


#'# Visualizing a genomic region with coMET




#' End of script