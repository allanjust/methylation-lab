#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  
#+ setdir02, echo = F
#' Local library
.libPaths("C:/EBC3/Rpackages")
knitr::opts_knit$set(root.dir = "../")

#' we have a processed dataset with 15 samples (otherwise we run script 01)
#if(!exists("WB.noob")){
#  source("code/meth01_process_data.R")
#}

# load the data
library(minfi)
load("C:/EBC3/Data/WB.noob.RData") # phenotype data
dim(WB.noob)
cellprop<-read.csv("C:/EBC3/methylation-lab/data/cellprop_WB_20samps_EPICdemo.csv") # cell type composition
load("C:/EBC3/Data/betas.rcp.RData") # processed betas
load("C:/EBC3/Data/Gbeta.RData") # annotation file

#' load packages
suppressPackageStartupMessages({
  library(CpGassoc) # for running association analysis between methylation levels values and phenotype of interest
  library(data.table) # for fast aggregation of large data 
  library(qqman) # for visualization of data
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) # for annotation for Illumina's EPIC methylation arrays
  library(DMRcate) # for regional analysis
  library(MASS) # for basic statistics
  library(sandwich) # for linear regression (robust sandwich variance estimator)
  library(lmtest) # for testing Linear Regression Models
  library(stringi) # string manipulation
})

#' consolidate our phenodata
pheno = data.frame(pData(WB.noob))
pheno = pheno[,c("SMOKE_STATUS","SEX","AGE","Array")]
pheno = cbind(pheno,cellprop)

pheno$SMOKE_STATUS = factor(pheno$SMOKE_STATUS)
pheno$SEX          = factor(pheno$SEX)

# cleanup 
rm(WB.noob)

#' quick check of the distribution of smoke between arrays
table(pheno[,c("SMOKE_STATUS","Array")])

## Cleaning up the methylation data
#' Filters a matrix of beta values by distance to single nucleotide polymorphism (SNP) and SNPs with the minor allele frequency (MAF) of 5% (rare variant). 
#' Also removes crosshybridising probes and sex-chromosome probes.
dim(betas.rcp)
betas.clean = rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG = dim(betas.clean)[1]
nCpG

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on Smoking status (as a predictor of methylation)  

#' First we can run a linear regression on a single CpG that we have already picked
CpG.name = "cg05575921"
CpG.level <- betas.clean[CpG.name,]

#' difference in methylation between smokers and non-smokers for this CpG
#' some descriptive statistics
knitr::kable(cbind(Min   = round( tapply(CpG.level,pheno$SMOKE_STATUS,min   ),3),
                   Mean  = round( tapply(CpG.level,pheno$SMOKE_STATUS,mean  ),3), 
                   Median= round( tapply(CpG.level,pheno$SMOKE_STATUS,median),3),
                   Max   = round( tapply(CpG.level,pheno$SMOKE_STATUS,max   ),3),
                   SD    = round( tapply(CpG.level,pheno$SMOKE_STATUS,sd    ),3),
                   N     = table( pheno$SMOKE_STATUS )))

#' difference in beta methylation values between Smokers and non smokers
par(mfrow=c(1,2))
boxplot(CpG.level ~ pheno$SMOKE_STATUS, main=paste0("Beta-values\n", CpG.name), col=c("blue","red"),ylim=c(.5,1))

#' linear regression on betas
summary(lm(CpG.level~pheno$SMOKE_STATUS))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' what if we use raw beta?
CpG.level.raw <- minfi::getBeta(WB)[CpG.name,] # select the same CpG

#' descriptive stats are different 
knitr::kable(cbind(Min    = round( tapply(CpG.level.raw, pheno$SMOKE_STATUS,min   ),3),
                   Mean   = round( tapply(CpG.level.raw, pheno$SMOKE_STATUS,mean  ),3), 
                   Median = round( tapply(CpG.level.raw, pheno$SMOKE_STATUS,median),3),
                   Max    = round( tapply(CpG.level.raw, pheno$SMOKE_STATUS,max   ),3),
                   SD     = round( tapply(CpG.level.raw, pheno$SMOKE_STATUS,sd    ),3),
                   N      = table( pheno$SMOKE_STATUS )))

boxplot(CpG.level.raw ~ pheno$SMOKE_STATUS, main=paste0("Raw Beta-values\n",CpG.name), col=c("blue","red"),ylim=c(.5,1))

#' linear regression on raw betas
summary(lm(CpG.level.raw ~ pheno$SMOKE_STATUS))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")] # results are different

#' comparison with m-values
CpG.mlevel = log2(CpG.level/(1-CpG.level))

knitr::kable(cbind(Min    = round( tapply(CpG.mlevel, pheno$SMOKE_STATUS,min   ),3),
                   Mean   = round( tapply(CpG.mlevel, pheno$SMOKE_STATUS,mean  ),3), 
                   Median = round( tapply(CpG.mlevel, pheno$SMOKE_STATUS,median),3),
                   Max    = round( tapply(CpG.mlevel, pheno$SMOKE_STATUS,max   ),3),
                   SD     = round( tapply(CpG.mlevel, pheno$SMOKE_STATUS,sd    ),3),
                   N      = table(pheno$SMOKE_STATUS)))

par(mfrow=c(1,2))
boxplot(CpG.level  ~ pheno$SMOKE_STATUS, main=paste0("Beta-values\n",CpG.name), col=c("blue","red"))
boxplot(CpG.mlevel ~ pheno$SMOKE_STATUS, main=paste0("M-values\n"   ,CpG.name), col=c("blue","red"))

#' linear regression on m-values
summary(lm(CpG.mlevel~pheno$SMOKE_STATUS))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' we can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
summary(lm(CpG.level~pheno$SMOKE_STATUS))$adj.r.squared

#' model on mvalues
summary(lm(CpG.mlevel~pheno$SMOKE_STATUS))$adj.r.squared

#'## EWAS and results using CpGassoc
#'see [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  

#' Smoking as predictor  
#' note that CpGassoc is quite fast for running almost a million regressions!

pheno$SMOKE = ifelse(pheno$SMOKE_STATUS=="SMOKER",1,0)
system.time(results1 <- cpg.assoc(betas.clean, pheno$SMOKE))

#' there are several components of the results
class(results1)
names(results1)
#' look at a few results  
#' here effect size is ~ mean difference in methylation proportion
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3]))
#' and the top hits
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3])[order(results1$results[,3]),])
#' check with previous result on our selected CpG (running lm without CpGassoc)
cbind(results1$coefficients[,4:5],results1$results[,c(1,3)])[CpG.name,]
summary(lm(CpG.level~pheno$SMOKE_STATUS))

#' Bonferroni significant hits
table(results1$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results1$results[,5] < 0.05)

#' EWAS with adjustment for cell types
#' now we can run the linear regression on betas adjusting for cell proportions
#' Need sex as indicator in covariate matrix
#' 
#' 
#' 
results2 = cpg.assoc(
           betas.clean
          ,pheno$SMOKE
          ,covariates=pheno[,c("SEX","AGE","CD8T","CD4T","NK","Bcell","Mono","Gran")]
          )

print(results2)

#'using mvalues
results3 <- cpg.assoc(
           betas.clean
          ,pheno$SMOKE
          ,covariates=pheno[,c("SEX","AGE","CD8T","CD4T","NK","Bcell","Mono","Gran")]
          ,logit.transform=TRUE
          )

print(results3)

#' ## Genomic inflation in EWAS
#' qqplot and lambda interpretation  
#+ fig.width=13, fig.height=7, dpi=300
par(mfrow=c(1,1))
plot(results1, main="QQ plot for association between methylation and Smoking\nadjusted for cell proportions")
plot(results2, main="QQ plot for association between (mvals) methylation and Smoking\nadjusted for cell proportions")

#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline
#' estimated at -log10(median=0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

#' Lambda before cell type adjustment
lambda(results1$results[,3])
#' Lambda after cell type adjustment
lambda(results2$results[,3])

#' Map the results to the epigenetic annotation
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))

#' Restrict to good quality probes and order data frames
IlluminaAnnot <- IlluminaAnnot[IlluminaAnnot$Name %in% results2$results$CPG.Labels,]
IlluminaAnnot <- IlluminaAnnot[match(results2$results$CPG.Labels, IlluminaAnnot$Name),]

#' Check that CpGs are align
identical(IlluminaAnnot$Name,results2$results$CPG.Labels)

datamanhat <- data.frame(CpG=results2$results[,1],Chr=as.character(IlluminaAnnot$chr),
                         Mapinfo=IlluminaAnnot$pos, UCSC_RefGene_Name=IlluminaAnnot$UCSC_RefGene_Name, 
                         Pval=results2$results[,3], Eff.Size = results2$coefficients[,4], Std.Error = results2$coefficients[,5])

#' see where the top hits are
head(datamanhat[order(datamanhat$Pval), ],n=7)

#' Volcano Plot-results2
#' with Bonferroni threshold and current FDR
plot(results2$coefficients[,4],-log10(results2$results[,3]), 
     xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot\nadjusted for cell proportions",ylim=c(0,8))
#Bonferroni threshold & FDR threshold
abline(h = -log10(0.05/(nCpG)), lty=1, col="#FDE725FF", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
#' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))

#' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
manhattan(datamanhat,"Chr","Mapinfo", "Pval", "CpG", 
          genomewideline = -log10(0.05/(nCpG)), suggestiveline = FALSE,
          main = "Manhattan Plot \n adjusted for cell proportions",ylim=c(0,8))

#' cleanup
rm(nCpG, CpG.name, CpG.level, CpG.mlevel, datamanhat, IlluminaAnnot,lambda,results1,results2,results3,Gbeta,WB,betas.rcp,CpG.level.raw);gc()
#' End of script 02
