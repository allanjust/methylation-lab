#'# Example of using a GEO dataset (450k)

library(utils)
library(data.table)
library(stringi)
library(minfi)

#+ setdir05, echo = F
knitr::opts_knit$set(root.dir = "../")

#' Investigating a 450k dataset from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) by NCBI
#' another great option is [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/browse.html) from the European EMBL-EBI 
#' 
#' Read about the dataset at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72556
#' more info here from the publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223358/

dir.create("data/GSE72556")

#' Either download manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72556 (look for GSE72556_RAW.tar under supplementary files)
#' or run this command in R
download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72556/suppl/GSE72556_RAW.tar",
              destfile="data/GSE72556/GSE72556_RAW.tar")

#' Extract the archive, again either manually, or by running 
untar("data/GSE72556/GSE72556_RAW.tar", exdir = "data/GSE72556")

#' Now download the meta data (list of samples and characteristics) from ...
#' https://gbnci-abcc.ncifcrf.gov/geo/
#' ... by searching for GSE72556 and retrieving the list of samples associated with this GEO series.
#' When you have the list of samples, click "Display options" and add "Characteristics Ch1" to the list of selected fields.
#' Download the .csv file and place it in the same folder.

#' Import the .csv file into R. The filename includes the current date, so edit the code below accordingly.
if(length(list.files(path = "data/GSE72556", pattern = "GEOmetadb_download.*\\.csv")) == 1){
    geometafile <- list.files(path = "data/GSE72556", pattern = "GEOmetadb_download.*\\.csv", full.names = TRUE)
  } else {
    warning(paste("did you download the metadata and put it in this directory?", file.path(getwd(), "data/GSE72556/")))
  }
meta <- fread(geometafile, sep=",")

if(!any(grepl("characteristics", names(meta)))){
  warning("characteristics field is missing - did it get added before download from GEOmetadb?")}

#' Parse the sample characteristics using regular expressions
meta[,adult_age   :=stri_match(characteristics_ch1,regex="adult age: (\\d+);"       )[,2]]
meta[,adult_bmi   :=stri_match(characteristics_ch1,regex="adult bmi: ([0-9\\.]+);"  )[,2]]
meta[,adult_waist :=stri_match(characteristics_ch1,regex="adult waist: ([0-9\\.]+);")[,2]]
meta[,child_gender:=stri_match(characteristics_ch1,regex="child gender: ([FM]);"    )[,2]]
meta[,child_age   :=stri_match(characteristics_ch1,regex="child age: (\\d);"        )[,2]]
meta[,child_bmi   :=stri_match(characteristics_ch1,regex="child bmi: ([0-9\\.]+);"  )[,2]]
meta[,child_waist :=stri_match(characteristics_ch1,regex="child waist: ([0-9\\.]+)$")[,2]]

#' Cast variables from type 'character' to 'numeric'
meta[,adult_age  :=as.numeric(adult_age)  ]
meta[,adult_bmi  :=as.numeric(adult_bmi)  ]
meta[,adult_waist:=as.numeric(adult_waist)]
meta[,child_age  :=as.numeric(child_age)  ]
meta[,child_bmi  :=as.numeric(child_bmi)  ]
meta[,child_waist:=as.numeric(child_waist)]

#' Determine the basenames of the .idat files required for minfi 'read.metharray'
meta[,basename:=paste0("data/GSE72556/", geo_accession, "_", title)]

#' Drop all the irrelevant columns
meta = meta[,list(adult_age,adult_bmi,adult_waist,child_gender,child_age,child_bmi,child_waist,basename)]

#' Drop rows with missing
sgithubmeta = na.omit(meta)

#' Import the methylation data
rgset = read.metharray(meta$basename)

#' Usually we should normalize the data, but we skip this step for now to save time
beta = getBeta(preprocessRaw(rgset)) 

#' Run a linear regression for the first 100 probes, store the p-values for 'adult_bmi'
pvals = apply(beta[1:100,],1,function(x){
	coef(summary(lm(x ~ adult_bmi+adult_age+adult_waist+child_gender+child_age+child_bmi+child_waist,data=meta)))['adult_bmi','Pr(>|t|)']
})

pvals
