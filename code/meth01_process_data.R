#' ---
#' title: "Processing, analyzing, and interpreting epigenome-wide DNA methylation data"
#' author: "Andrea Baccarelli, Andres Cardena, Elena Colicino, Cavin Ward-Caviness, Allan Just"
#' date: "June, 2017"
#' geometry: margin=2cm
#' number_sections: true
#' ---
#+ setdir01, echo = F
knitr::opts_knit$set(root.dir = "../")


#'# Processing methylation data  

#' Load packages that we will use focusing on *minfi*:  
#'  see [Aryee et al. Bioinformatics 2014](http://doi.org/10.1093/bioinformatics/btu049).  
#' Other popular options for processing & analyzing methylation data include RnBeads and methylumi
suppressMessages(library(minfi)) # popular package for methylation data
library(shinyMethyl) # for visualizing quality control
library(pryr) # for monitoring memory use
suppressMessages(library(matrixStats)) # for calculating summary statistics
library(ENmix) # probe type adjustment "rcp"
suppressMessages(library(limma)) # for MDS plots
suppressMessages(library(reshape,scales)) # reshape data and graphig 
suppressMessages(require(sva)) # for addressing batch effects
library(IlluminaHumanMethylationEPICmanifest) # manifest for Illumina's EPIC methylation 
library(IlluminaHumanMethylationEPICanno.ilmn10b2.hg19) # annotation for Illumina's EPIC methylation arrays.
#' Here we will use an existing methylation dataset - EPIC in peripheral blood (20 samples)  
#' To read in your own dataset (usually "idat files" ending in .idat)  
#' See the help for ?read.metharray.exp  
#library(EPICdemo) # example dataset - not a publically available package

#' import EPIC data from a sample sheet and idat files
idatPath <- "~/BootCamp_Epigenetics/Data/EPICdemo/inst/extdata"
targets <- read.csv("~/BootCamp_Epigenetics/Data/EPICdemo/inst/extdata/sample.info.csv", strip.white=T, stringsAsFactors=F)
targets$Basename <- paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)
WB <- read.metharray.exp(base=idatPath, targets=targets, verbose=T)
ncol(WB)
#' alternative way to import EPIC data, using EPICdemo
#sheet <- read.metharray.sheet(base = system.file("extdata", package = "EPICdemo"), pattern = "csv$")
#WB <- read.metharray.exp(targets = sheet)


#' Look at the attributes of this dataset  
#' It is stored as a RGChannelSet which means it is not yet processed (red and green signals stored separately)
WB

#' Prepare a summary of the dataset to look at interactively
summaryqc <- shinySummarize(WB)
#' Visualize qc with shinyMethyl
#runShinyMethyl(summaryqc)
#' cleanup
#rm(sheet, covariates, colorSet, current.control.type, current.density.type, current.probe.type, 
#   first_time, genderCutoff, mouse.click.indices, sampleColors, summaryqc)

#' Look at some of the phenotype data:
pData(WB)[,1:7]

#' ## Estimate cell proportions
#' Estimating proportions of 6 cell types found in peripheral blood 
#'  _This next command is commented out because it requires >4GB of RAM_  
#'  if you don't have that - you can load the presaved output below
# cellprop <- estimateCellCounts(WB, compositeCellType = "Blood",
#   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
# write.csv(cellprop, file = "data/cellprop_WB_20samps_EPICdemo.csv", row.names = F)
#' read in the estimated cell proportions:
cellprop <- read.csv("data/cellprop_WB_20samps_EPICdemo.csv")
#' Here are the estimates
#+ cellprop, results = "asis"
knitr::kable(cellprop, digits = 2)
#' note that they are close to summing to 1
summary(rowSums(cellprop))
#'Distribution of estimated cell types
boxplot(cellprop*100, col=1:ncol(cellprop),xlab="Cell type",ylab="Estimated %",main="Cell type distribution")
#'Distribution of estimated cell types by smoking status
meltData <- melt(cbind(cellprop*100,pData(WB)$SMOKE_STATUS))
names(meltData)[1:3]<-c('Smoking','Celltype','Proportion')
boxplot(Proportion ~Smoking+Celltype, data=meltData,col=c("blue","red"),xaxt="n",main="Cell type distribution by smoking status")
axis(1,at=seq(from=1.5, to=11.5,by=2),adj=1,labels=c(colnames(cellprop)))
legend("topleft",c("Non-smoker","Smoker"),pch=15,bty='n',col=c("blue","red"),cex=1.3)

#'
#' ## Preprocessing the data
#' Preprocess the data - this removes technical variation  
#' There are several popular methods including intra- and inter-sample normalizations  
#' Here we demonstrate an effective and lightweight approach:  
#' "Normal out of band background" (Noob) within-sample correction - see [Triche et al 2013](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627582/)
system.time(WB.noob <- preprocessNoob(WB))
#' We see the resulting object is now a MethylSet (because the RGset has been preprocessed)  
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
densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
densityPlot(WB.noob, add = F, pal = "#FDE725FF")
# Add legend
legend("topleft", c("Noob","Raw"), 
       lty=c(1,1), title="Normalization", 
       bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))
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
#' Let's look at the median detection P-values
#+ fig.width=8, fig.height=6, dpi=300
barplot(colMeans(detect.p), col=rainbow(dim(detect.p)[2]), las=2, 
        cex.names=0.7, main="Mean detection P by sample",cex.axis=0.8, ylim=c(0,3e-4))
#' Count of failed probes per sample 
#' P>0.01 (i.e. not significant compared to background signal)
knitr::kable(t(as.matrix(colSums(detect.p > 0.01))[1:5,]))
knitr::kable(t(as.matrix(colSums(detect.p > 0.01))[6:10,]))
knitr::kable(t(as.matrix(colSums(detect.p > 0.01))[11:15,]))
knitr::kable(t(as.matrix(colSums(detect.p > 0.01))[16:20,]))
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
rm(intersect, detect.p)

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
# almost 84% of our probes are type II
knitr::kable(t(table(onetwo)))

#' Density plots by Infinium type: before and after RCP calibration  
#' Probe-type bias adjustment before and after RCP
#+ fig.width=15, fig.height=7, dpi=300
par(mfrow=c(1,2)) # Side-by-side density distributions 
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density',ylim=c(0,6.5))
densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density probe-type adjusted',ylim=c(0,6.5))
densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
legend("center", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n',col=c("#FDE725FF","#440154FF"))
#' notice that the type I and II peaks are more closely aligned after rcp adjustment  
#' (particularly in the higher peak)
rm(onetwo, typeI, typeII)

#' ## Batch effects
#' As an example of an observable batch effect, samples are processed in plates (e.g. bisulfite converting 96 at a time),  
#' and then in chips (EPIC array has 8 rows and 1 column).  
#' This can create batch effects (technical variation) with different intensities by position (row effect).  
#' Other commonly observed batch effects include bisulfite processing plate, chip, and processing date.
#' Let's check if samples varied across rows in these data (even in different chips):
knitr::kable(t(as.matrix(table(pData(WB.noob)$Array))), 
             col.names = paste0("Row ", 1:8))

#' ## Principal Component Analysis (PCA)
#' Calculate major sources of variability of DNA methylation using PCA
PCobject = prcomp(t(betas.rcp), retx = T, center = T, scale. = T)

#' Extract the Principal Components from SVD
PCs <- PCobject$x
#' Proportion of variance explained by each additional PC
cummvar <- summary(PCobject)$importance["Cumulative Proportion", 1:10]
knitr::kable(t(as.matrix(cummvar)),digits = 2)

#' Is the major source of variability associated with position on chip?
par(mfrow = c(1, 1))
boxplot(PCs[, 1] ~ pData(WB.noob)$Array,
        ylab = "PC1",las=2, main="Position on chip")
summary(lm(PCs[, 1] ~ pData(WB.noob)$Array))

#' ## Removing batch effects using ComBat from the sva package
# First we convert from beta-values to M-values
Mvals <- log2(betas.rcp)-log2(1-betas.rcp)
#' ComBat eBayes adjustment using a known variable of interest (here we use row)
Mvals.ComBat <- ComBat.mc(Mvals, batch = pData(WB.noob)$Array,nCores = detectCores()-1)
# Convert M-values back to beta-values
betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)

#' PCA after removing batch effects
PCobject <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
PCs <- PCobject$x
#' The first PC is no longer associated with sample plate
boxplot(PCs[, 1] ~ pData(WB.noob)$Array,
        ylab = "PC1",las=2, main="Position on chip")
summary(lm(PCs[, 1] ~ pData(WB.noob)$Array))
#' ComBat removed the apparent batch effect
#cleanup
rm(PCs, Mvals, cummvar, PCobject)

#' Phenotype associated with large sources of variability
#+ fig.width=8, fig.height=6, dpi=300
plotMDS(betas.rcp, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pData(WB.noob)$SEX)],
        dim=c(1,2),cex=1.5)
legend("topright", legend=levels(factor(pData(WB.noob)$SEX)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))


#'## predict sex from methylation
Gbeta <- mapToGenome(WB)  #map to the genome
#' getSex predicts sex based on X and Y chromosome methylation intensity
getSex(Gbeta) 
#' we see that our predictions match the phenodata
table(pData(WB)$SEX,getSex(Gbeta)$predictedSex)
#' We can actually look at the intensities in the sex-crhomosomes
plotSex(getSex(Gbeta, cutoff = -2),id=pData(WB)$SEX)
# cleanup


#' # memory usage
#' as a reminder - these are large datasets and we are working in RAM.  
#' Check memory usage with:
pryr::mem_used()
# this command will check the maximum memory usage on Windows
memory.size(max = T)

#' End of script 1
#' 
