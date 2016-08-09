#' ---
#' title: "Processing, analyzing, and interpreting epigenome-wide DNA methylation data"
#' author: "Andrea Baccarelli, Andres Cardena, Elena Colicino, Allan Just"
#' date: "September, 2016"
#' geometry: margin=2cm
#' number_sections: true
#' ---

#'# Processing methylation data  
#'
#' Here we will use an existing methylation dataset - 450k in cord blood (15 samples)  
#' To read in your own dataset (usually "idat files" ending in .idat)  
#' See the help for ?read.metharray.exp  

#'
#' Load packages that we will use focusing on *minfi*:  
#'  see Aryee et al. Bioinformatics 2014. http://doi.org/10.1093/bioinformatics/btu049  
#' Other popular options for processing & analyzing methylation data include RnBeads and methylumi
suppressMessages(library(minfi)) # popular package for methylation data
library(FlowSorted.CordBlood.450k) # example dataset
library(pryr) # for monitoring memory use
library(ENmix) # probe type adjustment "rcp"

#' bring a 450k dataset in to memory (from eponymous BioC package above)  
#'   see Bakulski et al. Epigenetics 2016.
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
cellprop <- estimateCellCounts(WB, compositeCellType = "CordBlood",
  cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
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
#' "Normal out of band background" (Noob) within-sample correction - see Triche et al 2013
system.time(WB.noob <- preprocessNoob(WB))
#' We see the resulting object is now a MethylSet (because the RGset has been preprocessed)
WB.noob


#' we can look at a few methylation values on the fly
#' and see that the preprocessing changed them
# first three CpGs on the first three samples
# raw RGset
print(getBeta(WB)[1:3,1:3], digits = 2)
# after preprocessing
print(getBeta(WB.noob)[1:3,1:3], digits = 2)

#' Distribution of beta-values: before and after normalization
#+ fig.width=8, fig.height=6, dpi=300
densityPlot(WB, main = "density plots before and after preprocessing", pal="blue")
densityPlot(WB.noob, add = F, pal = "magenta")
# Add legend
legend("topright", c("Noob","Raw"), 
lty=c(1,1), title="Normalization", 
bty='n', cex=0.8, col=c("blue","magenta"))

#' ### probe failures due to low intensities
#' We want to drop probes that were not significantly above background signal (from negative control probes)
detect.p <- detectionP(WB, type = "m+u")
#' Proportion of failed probes per sample 
#' P>0.01 (i.e. not significant compared to background signal)
knitr::kable(t(as.matrix(colMeans(detect.p > 0.01) * 100)), digits = 2)
#' Number of failed probes
sum(detect.p > 0.01)
#'Restrict data to good probes only
detect.p[detect.p > 0.01] <- NA
detect.p <- na.omit(detect.p)
intersect <- intersect(rownames(getAnnotation(WB)), rownames(detect.p))
length(intersect)
#' Filter bad probes from our methylset
length(rownames(WB.noob))
WB.noob <- WB.noob[rownames(getAnnotation(WB.noob)) %in% intersect,]
length(rownames(WB.noob))
# cleanup
rm(intersect, detect.p, WB)

#' Need to adjust for probe-type bias Infinium I (type I) and Infinium II (type II) probes
## RCP with EnMix: Regression on Correlated Probes (Niu et al. Bioinformatics 2016)
betas.rcp <- rcp(WB.noob)
dim(betas.rcp)
#' note that this package takes beta values out of the minfi object - result is a matrix
class(betas.rcp) 

## Annotation of Infinium type for each probe (I vs II)
typeI <-   minfi::getProbeInfo(WB.noob,type="I")$Name
typeII <-  minfi::getProbeInfo(WB.noob,type="II")$Name
onetwo <- rep(1, nrow(betas.rcp))
onetwo[rownames(betas.rcp) %in% typeII] <- 2
knitr::kable(t(table(onetwo)))

#' Density plots by Infinium type: before and after RCP calibration
#+ fig.width=14, fig.height=7, dpi=300
par(mfrow=c(1,2)) # Side-by-side density distributions 
# Noob adjusted Beta density
plot(density(getBeta(WB.noob)[,1][onetwo==1]), col="blue",ylim=c(0,6), 
     main='Beta density',xlab=expression(beta~"-value")) # plot the first density
for(i in 2:dim(getBeta(WB.noob))[2]){          # Add the lines to the existing plot
  lines(density(getBeta(WB.noob)[,i][onetwo==1]), col="blue")        
}
for(i in 1:dim(getBeta(WB.noob))[2]){          # Add the lines to the existing plot
  lines(density(getBeta(WB.noob)[,i][onetwo==2]), col="red")        
}
legend("topright", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n', cex=0.8, col=c("blue","red"))

#' RCP adjusted Beta Density
plot(density(betas.rcp[,1][onetwo==1]), col="blue",ylim=c(0,6), 
     main='Beta density probe-type adjusted',xlab=expression(beta~"-value")) # plot the first density
for(i in 2:dim(betas.rcp)[2]){          # Add the lines to the existing plot
  lines(density(betas.rcp[,i][onetwo==1]), col="blue")        
}
for(i in 1:dim(betas.rcp)[2]){          # Add the lines to the existing plot
  lines(density(betas.rcp[,i][onetwo==2]), col="red")        
}
legend("topright", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n', cex=0.8, col=c("blue","red"))
rm(i, onetwo, typeI, typeII)

#' ## Batch effects
#' Samples are processed in plates  
#' this can create batch effects with different intensities by plate  
#' Let's check if we have a plate effects in our data:
knitr::kable(t(as.matrix(table(pData(WB.noob)$Plate_ID))), col.names = c("Plate 1","Plate2"))


#' ## Principal Component Analysis for the DNA methylation data
#' Calculate major sources of variability of DNA methylation using PCA
PCobject = prcomp(t(betas.rcp), retx = T, center = T, scale. = T)

#' Extract the Principal Components from SVD
PCs <- PCobject$x
#' Proportion of variance explained by each PC
cummvar <- summary(PCobject)$importance["Cumulative Proportion", 1:10]
knitr::kable(t(as.matrix(cummvar)),digits = 2)


#' Is the major source of variability associated with sample plate?
boxplot(PCs[,1]~pData(WB.noob)$Plate_ID,
        xlab = "Sample Plate",ylab="PC1",
        col=c("red","blue"))
t.test(PCs[,1]~pData(WB.noob)$Plate_ID)

#' Removing batch effects using ComBat from the sva package
suppressMessages(require(sva))
# From Beta-values to M-values
Mvals <- log2(betas.rcp)-log2(1-betas.rcp)
Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$Plate_ID)
# From M-values back to Beta-values
betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)

#' PCA after removing batch effects
PCobject <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
PCs <- PCobject$x
cummvar <- summary(PCobject)$importance["Cumulative Proportion", 1:10]
knitr::kable(t(as.matrix(cummvar)),digits = 2)
#' The first PC is no longer associated with sample plate
boxplot(PCs[,1] ~ pData(WB.noob)$Plate_ID,
        xlab = "Sample Plate", ylab = "PC1",
        col = c("red","blue"))
t.test(PCs[,1] ~ pData(WB.noob)$Plate_ID)
#' ComBat removed the apparent batch effect
#cleanup
rm(PCs, Mvals, cummvar, PCobject)

#' just checking how much memory we are using
#mem_used()

#' this command will check the maximum memory usage on Windows
memory.size(max = T)

#' End of script 1
#' 