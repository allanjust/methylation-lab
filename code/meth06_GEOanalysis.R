#'# Example of using a GEO dataset (450k)

.libPaths("C:/EBC4/Rpackages")

options(warn=-1)
library(utils)
library(stringi)
suppressMessages(library(minfi))
suppressMessages(library(purrr))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
options(warn=0)

#' Investigating a 450k dataset from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) by NCBI
#' another great option is [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/browse.html) from the European EMBL-EBI 
#' 
#' Read about the dataset at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72556
#' more info here from the publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223358/

dir.create("GSE72556")

#' Download meta data (copy URL in browser to see what the requested file looks like)
meta = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72556&targ=gsm&form=text&view=brief"
meta = readLines(meta)

#' Split into individual samples
meta = split(meta,cumsum(meta %like% "^\\^SAMPLE = GSM"))

#' Extract GSM accessions
names(meta) = map(meta,1) %>% stri_match_first(regex="GSM\\d+")

#' Parse meta data
imap(meta,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> meta

meta = rbindlist(meta)

#' Keep only information on sample characteristics and supplementary files
meta = meta[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = meta[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = meta$value[i] %>% stri_split(fixed=": ")
meta$variable[i] = map_chr(ch,1)
meta$value   [i] = map_chr(ch,2)
rm(ch)

#' Find the URLs pointing to the two .idat files
meta[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
meta[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

#' Reshape data.table from long to wide format
meta = dcast(meta, gsm ~ variable)

#' Keep only the first ten samples
meta = meta[1:10]

#' Download .idat files
map2(meta$red, "GSE72556/" %s+% meta$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(meta$grn, "GSE72556/" %s+% meta$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
meta$red = NULL; meta$grn = NULL

#' Type casting
meta$`adult age`    %<>% as.numeric
meta$`adult bmi`    %<>% as.numeric
meta$`adult waist`  %<>% as.numeric
meta$`child age`    %<>% as.numeric
meta$`child bmi`    %<>% as.numeric
meta$`child waist`  %<>% as.numeric
meta$`child gender` %<>% factor

#' Import the methylation data
rgset = read.metharray("GSE72556/" %s+% meta$gsm)

#' Usually we should normalize the data, but we skip this step for now to save time
beta = getBeta(preprocessRaw(rgset)) 

# Run a linear regression for the first 100 probes, store the p-values for 'adult_bmi'
pvals = apply(beta[1:100,],1,function(x){
	coef(summary(lm(x ~ `adult bmi`+`adult age`+`adult waist`+`child gender`+`child age`+`child bmi`+`child waist`,data=meta)))["`adult bmi`","Pr(>|t|)"]
})

pvals
