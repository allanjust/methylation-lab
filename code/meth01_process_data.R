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
#' See our additional script "meth00_import_idat.R" *not yet written*

#'
#' Load packages that we will use focusing on *minfi*:  
#'  see Aryee et al. Bioinformatics 2014. http://doi.org/10.1093/bioinformatics/btu049
#' Other popular options for processing & analyzing methylation data include RnBeads and methylumi
suppressMessages(library(minfi)) # popular package for methylation data
library(FlowSorted.CordBlood.450k) # example dataset
library(pryr) # for monitoring memory use

#' bring a 450k dataset in to memory (from eponymous BioC package above)
#'   see Bakulski et al. Epigenetics 2016.
data(FlowSorted.CordBlood.450k)

#' because it is flow sorted - the authors give us the cell types  
#' here we show the frequencies
table(pData(FlowSorted.CordBlood.450k)$CellType)

# subset to just the Whole Blood samples since this is the most common for epi studies
WB <- FlowSorted.CordBlood.450k[, pData(FlowSorted.CordBlood.450k)$CellType == "WholeBlood"]
ncol(WB)
#' drop the reference dataset from memory
rm(FlowSorted.CordBlood.450k)
#' we need to change the sampleName attribute here just because we are using a reference and
#' these samples are otherwise in the reference dataset when we want to estimate their composition.
#' #' your personal samples won't need to be renamed
sampleNames(WB) <- 1:15
#' Look at the attributes of this dataset
#' It is stored as a RGChannelSet which means it is not yet processed (red and green signals stored separately)
WB

#'
#' ## Estimate cell proportions
#' Estimating proportions of 7 cell types found in cord blood (note the nucleated Red Blood Cells)
cellprop <- estimateCellCounts(WB, compositeCellType = "CordBlood",
  cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
#' Here are the estimates
#+ cellprop, results = "asis"
knitr::kable(cellprop, digits = 2)
#' note that they are close to summing to 1
summary(rowSums(cellprop))
#'Distribution of estiamted cell types
boxplot(cellprop*100, col=1:ncol(cellprop),xlab="Cell type",ylab="Estimated %")

#'
#' ## Preprocessing the data
#' Preprocess the data - this removes technical variation  
#' There are several popular methods including intra- and inter-sample normalizations  
#' Here we demonstrate an effective and lightweight approach
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
plot(density(getBeta(WB.noob)[,1]), col="blue",ylim=c(0,5),
     main='Beta Density',xlab=expression(beta~"-value")) # plot the first density
for(i in 2:dim(getBeta(WB.noob))[2]){ #Plot remaning noob adjusted samples
  lines(density(getBeta(WB.noob)[,i]), col="blue")      
}
for(i in 1:dim(getBeta(WB.noob))[2]){ # Add lines to the existing plot for unadjusted betas
  lines(density(na.omit(getBeta(WB))[,i]), col="magenta")        
}
legend("topright", c("Noob","Raw"), 
       lty=c(1,1), title="Normalization", 
       bty='n', cex=0.8, col=c("blue","magenta"))

#' Need to adjust for probe-type bias Infinium I (type I) and Infinium II (type II) probes
## BMIQ with EnMix: multi-processor wrapper of BMIQ
require(ENmix)
cores<-detectCores()
betas.bmiq<-bmiq.mc(WB.noob, nCores=cores-1) #4-cores~2.8 min;3-cores~3.3 min
dim(betas.bmiq)

## Annotation of Infinium type for each probe (I vs II)
typeI <-   minfi::getProbeInfo(WB.noob,type="I")$Name
typeII <-  minfi::getProbeInfo(WB.noob,type="II")$Name
onetwo <- rep(1, nrow(betas.bmiq))
onetwo[rownames(betas.bmiq) %in% typeII] <- 2
table(onetwo)


#' Density plots by Infinium type: before and after BMIQ 
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

# BMIQ adjusted Beta Density
plot(density(betas.bmiq[,1][onetwo==1]), col="blue",ylim=c(0,6), 
     main='Beta density BMIQ adjusted',xlab=expression(beta~"-value")) # plot the first density
for(i in 2:dim(betas.bmiq)[2]){          # Add the lines to the existing plot
  lines(density(betas.bmiq[,i][onetwo==1]), col="blue")        
}
for(i in 1:dim(betas.bmiq)[2]){          # Add the lines to the existing plot
  lines(density(betas.bmiq[,i][onetwo==2]), col="red")        
}
legend("topright", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n', cex=0.8, col=c("blue","red"))

#' just checking how much memory we are using
#mem_used()

#' this command will check the maximum memory usage on Windows
memory.size(max = T)

#' End of script