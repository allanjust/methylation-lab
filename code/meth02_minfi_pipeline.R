#'# minfi-ewastools pipeline comparison
.libPaths("C:/EBC4/Rpackages")
#+ setdir01, echo = F
knitr::opts_knit$set(root.dir = "../")

options(warn=-1)
library(data.table)
library(magrittr)
options(warn=0)

#' List of idat files
pheno = fread("C:/EBC4/methylation-lab/data/pheno_clean.csv")

#' ## Pre-processing with `ewastools`

library(ewastools)

beta =
	paste0("data/",pheno$gsm) %>% # idat file paths
	read_idats(quiet=TRUE) %>%    # import fluorescence intensities
	detectionP %>%                # compute detection p-values
	mask(0.01) %>%                # set undetected probes to missing
	correct_dye_bias %>%
	dont_normalize                # calculate beta-values

dim(beta)

detach("package:ewastools")

#' ## Pre-processing with `minfi`

options(warn=-1)
suppressMessages(library(minfi))
suppressMessages(library(ENmix))
options(warn=0)

idat_files = paste0("data/",pheno$gsm)

# importing idat files, result is a RGChannelSet
rgset = read.metharray(basenames=idat_files)
rgset

# `preprocessRaw` computes beta-value based on raw flourescence intensities
methylset = preprocessRaw(rgset)
methylset

#' Many other preprocessing methods are available
#' 
#' * `preprocessNoob`: background substraction/correction
#' * `preprocessFunnorm`: optionally includes preprocessNoob
#' * `preprocessIllumina`: normalization as in GenomeStudio software
#' * `preprocessQuantile`: quantile normalization stratified by U/M signal and probe type design
#' * `preprocessSWAN`: quantile normalization stratified by probe type design and number of CpG sites in probe sequence
#' * `preprocessENmix`: background correction (`ENmix` package)

enmix = preprocessENmix(rgset[,1:2])

#' Due to their design, Type II probes feature more background noise, shifting the peaks for completely (un)methylated Cpg sites towards 0.5
plotBetasByType(getBeta(enmix)[,1],probeTypes=getAnnotation(enmix))

#' Peaks of Type I and Type II beta-value distributions are now aligned
plotBetasByType(    rcp(enmix)[,1],probeTypes=getAnnotation(enmix))


#' ## Cell type prediction
#' 
minfi.LC = estimateCellCounts(rgset,compositeCellType="Blood")

ewastools.LC = ewastools::estimateLC(beta,ref="Reinius")

cor(ewastools.LC$CD4,minfi.LC[,"CD4T"] )
cor(ewastools.LC$GR ,minfi.LC[,"Gran"] )
cor(ewastools.LC$MO ,minfi.LC[,"Mono"] )
cor(ewastools.LC$B  ,minfi.LC[,"Bcell"])

#' In case you have samples measured on both the 450K and EPIC chips, you can virtually convert them
convertArray(rgset,outType="IlluminaHumanMethylationEPIC")


#' 'ewastools' implements a more stringent approach to calculate detection p-values
#' 
detP1 = detectionP(rgset[,5]) # Sample #5 is female
detP2 = ewastools::detectionP.minfi(rgset[,5])

chrY = getAnnotation(methylset)$chr == "chrY"

#' Most Y chromosome probes are called detected using the minfi approach
table(detP1[chrY,] < 0.05)

#' ... whereas `detectionP.minfi` from `ewastools` results in a more accurate classification
table(detP2[chrY,] < 0.05)
