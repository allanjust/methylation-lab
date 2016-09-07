#' ---
#' title: "Processing, analyzing, and interpreting epigenome-wide DNA methylation data"
#' author: "Andrea Baccarelli, Andres Cardena, Elena Colicino, Allan Just"
#' date: "September, 2016"
#' geometry: margin=2cm
#' number_sections: true
#' ---
#+ setdir01, echo = F
knitr::opts_knit$set(root.dir = "../")


#'# Processing methylation data  
#'
#' Here we will use an existing methylation dataset - 450k in cord blood (15 samples)  
#' To read in your own dataset (usually "idat files" ending in .idat)  
#' See the help for ?read.metharray.exp  

#'
#' Load packages that we will use focusing on *minfi*:  
#'  see [Aryee et al. Bioinformatics 2014](http://doi.org/10.1093/bioinformatics/btu049).  
#' Other popular options for processing & analyzing methylation data include RnBeads and methylumi
suppressMessages(library(minfi)) # popular package for methylation data
library(FlowSorted.CordBlood.450k) # example dataset
library(pryr) # for monitoring memory use
suppressMessages(library(matrixStats)) # for calculating summary statistics
library(ENmix) # probe type adjustment "rcp"
suppressMessages(require(sva)) # for addressing batch effects

#' bring a 450k dataset in to memory (from eponymous BioC package above)  
#'   see [Bakulski et al. Epigenetics 2016](http://www.ncbi.nlm.nih.gov/pubmed/27019159).
data(FlowSorted.CordBlood.450k)

#' because it is flow sorted - the authors give us the cell types  
#' here we show the frequencies:
table(pData(FlowSorted.CordBlood.450k)$CellType)

# subset to just the Whole Blood samples since this is the most common for epi studies
WB <- FlowSorted.CordBlood.450k[, pData(FlowSorted.CordBlood.450k)$CellType == "WholeBlood"]
ncol(WB)
#' we need to change the sampleName attribute here just because we are using a reference and
#' these samples are otherwise in the reference dataset when we want to estimate their composition.
# your personal samples won't need to be renamed
sampleNames(WB) <- 1:15
#' Look at the attributes of this dataset  
#' It is stored as a RGChannelSet which means it is not yet processed (red and green signals stored separately)
WB

#'
#' ## Estimate cell proportions
#' Estimating proportions of 7 cell types found in cord blood (note the nucleated Red Blood Cells)  
#'  _This next command is commented out because it requires >4GB of RAM_  
#'  if you don't have that - you can load the presaved output below
# cellprop <- estimateCellCounts(WB, compositeCellType = "CordBlood",
#   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
# write.csv(cellprop, file = "data/cellprop_WB_15samps_bakulski2016.csv", row.names = F)
#' read in the estimated cell proportions:
cellprop <- read.csv("data/cellprop_WB_15samps_bakulski2016.csv")
#' drop the reference dataset from memory
rm(FlowSorted.CordBlood.450k)
#' Here are the estimates
#+ cellprop, results = "asis"
knitr::kable(cellprop, digits = 2)
#' note that they are close to summing to 1
summary(rowSums(cellprop))
#'Distribution of estimated cell types
boxplot(cellprop*100, col=1:ncol(cellprop),xlab="Cell type",ylab="Estimated %")

#'
#' ## Preprocessing the data
#' Preprocess the data - this removes technical variation  
#' There are several popular methods including intra- and inter-sample normalizations  
#' Here we demonstrate an effective and lightweight approach:  
#' "Normal out of band background" (Noob) within-sample correction - see [Triche et al 2013](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627582/)
system.time(WB.noob <- preprocessNoob(WB))
#' We see the resulting object is now a MethylSet (because the RGset has been preprocessed)  
#' Minfi is incorrectly saying the data are still raw - we verify this is not true below
WB.noob


#' we can look at a few methylation values on the fly
#' and see that the preprocessing changed them:
# first three CpGs on the first three samples
# raw RGset
print(getBeta(WB)[1:3,1:3], digits = 2)
# after preprocessing
print(getBeta(WB.noob)[1:3,1:3], digits = 2)
#' and we can check that
#' the dimensions of the beta matrices are unchanged
identical(dim(getBeta(WB)), dim(getBeta(WB.noob)))

#' Distribution of beta-values: before and after normalization
#+ fig.width=8, fig.height=6, dpi=300
densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF")
densityPlot(WB.noob, add = F, pal = "#FDE725FF")
# Add legend
legend("topright", c("Noob","Raw"), 
  lty=c(1,1), title="Normalization", 
  bty='n', cex=0.8, col=c("#FDE725FF","#440154FF"))
#' notice the blue density traces (raw) are more spread out; background correction brings them together  

#' the use of a density plot may give the impression of values outside (0,1)  
#' but let's check:
print(colRanges(getBeta(WB), na.rm = T), digits = 2)
print(colRanges(getBeta(WB.noob)), digits = 2)
#' we confirm that plot tails were an artifact of applying a continuous distribution  
#'  (in estimating the density function)

#' ## probe failures due to low intensities  
#' We want to drop probes with intensity that is not significantly above background signal (from negative control probes)
detect.p <- detectionP(WB, type = "m+u")
#' Count of failed probes per sample 
#' P>0.01 (i.e. not significant compared to background signal)
knitr::kable(t(as.matrix(colSums(detect.p > 0.01))))
#' This is less than 1% of probes per sample  
#' Total number of failed measures (in all probes)
sum(detect.p > 0.01)
#'Restrict data to good probes only:
detect.p[detect.p > 0.01] <- NA
detect.p <- na.omit(detect.p)
intersect <- intersect(rownames(getAnnotation(WB)), rownames(detect.p))
length(intersect)
#' Filter bad probes from our methylset
nrow(WB.noob)
WB.noob <- WB.noob[rownames(getAnnotation(WB.noob)) %in% intersect,]
nrow(WB.noob)
# cleanup
rm(intersect, detect.p, WB)

#' #Probe type adjustment  
#' Need to adjust for probe-type bias Infinium I (type I) and Infinium II (type II) probes  
#' RCP with EnMix: Regression on Correlated Probes [Niu et al. Bioinformatics 2016](http://www.ncbi.nlm.nih.gov/pubmed/27153672)
betas.rcp <- rcp(WB.noob)
dim(betas.rcp)
#' note that this package takes beta values out of the minfi object - result is a matrix
class(betas.rcp) 

## Annotation of Infinium type for each probe (I vs II)
typeI <-   minfi::getProbeInfo(WB.noob,type="I")$Name
typeII <-  minfi::getProbeInfo(WB.noob,type="II")$Name
onetwo <- rep(1, nrow(betas.rcp))
onetwo[rownames(betas.rcp) %in% typeII] <- 2
# almost three quarter of our probes are type II
knitr::kable(t(table(onetwo)))

#' Density plots by Infinium type: before and after RCP calibration  
#' Probe-type bias adjustment before and after RCP
#+ fig.width=15, fig.height=7, dpi=300
par(mfrow=c(1,2)) # Side-by-side density distributions 
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density')
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density probe-type adjusted')
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
legend("topright", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n',col=c("#FDE725FF","#440154FF"))
#' notice that the type I and II peaks are more closely aligned after rcp adjustment  
#' (particularly in the higher peak)
rm(onetwo, typeI, typeII)

#' ## Batch effects
#' As an example of an observable batch effect, samples are processed in plates (e.g. bisulfite converting 96 at a time).  
#' This can create batch effects (technical variation) with different intensities by plate.  
#' Other commonly observed batch effects include the position on the chip (e.g. the row effect).  
#' Let's check if samples were on different plates in these data:
knitr::kable(t(as.matrix(table(pData(WB.noob)$Plate_ID))), col.names = c("Plate 1","Plate2"))


#' ## Principal Component Analysis (PCA)
#' Calculate major sources of variability of DNA methylation using PCA
PCobject = prcomp(t(betas.rcp), retx = T, center = T, scale. = T)

#' Extract the Principal Components from SVD
PCs <- PCobject$x
#' Proportion of variance explained by each additional PC
cummvar <- summary(PCobject)$importance["Cumulative Proportion", 1:10]
knitr::kable(t(as.matrix(cummvar)),digits = 2)


#' Is the major source of variability associated with sample plate?
par(mfrow = c(1, 1))
boxplot(PCs[, 1] ~ pData(WB.noob)$Plate_ID,
        xlab = "Sample Plate", ylab = "PC1",
        col = c("#FDE725FF", "#440154FF"))
t.test(PCs[, 1] ~ pData(WB.noob)$Plate_ID)

#' ## Removing batch effects using ComBat from the sva package
# First we convert from beta-values to M-values
Mvals <- log2(betas.rcp)-log2(1-betas.rcp)
#' ComBat eBayes adjustment using a known variable of interest (here we use plate)
Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$Plate_ID)
# Convert M-values back to beta-values
betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)

#' PCA after removing batch effects
PCobject <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
PCs <- PCobject$x
#' The first PC is no longer associated with sample plate
boxplot(PCs[,1] ~ pData(WB.noob)$Plate_ID,
        xlab = "Sample Plate", ylab = "PC1",
        col = c("#FDE725FF","#440154FF"))
t.test(PCs[,1] ~ pData(WB.noob)$Plate_ID)
#' ComBat removed the apparent batch effect
#cleanup
rm(PCs, Mvals, cummvar, PCobject)

#' # memory usage
#' as a reminder - these are large datasets and we are working in RAM.  
#' Check memory usage with:
pryr::mem_used()
# this command will check the maximum memory usage on Windows
memory.size(max = T)

#' End of script 1
#' 