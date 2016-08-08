#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  


#' we have a processed dataset with 15 samples (otherwise we run script 01)
if(!exists("WB.noob")){
  rmarkdown::render("meth01_process_data.R")
}
dim(WB.noob)

#libraries
suppressPackageStartupMessages({
  library(CpGassoc)
  library(data.table)
  library(qqman)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(DMRcate)
})

## predict gender from methylation
Gbeta<-mapToGenome(WB.noob)  #map to the genome
getSex(Gbeta) #get sex
cbind(Sex=pData(WB.noob)$Sex,PSex=getSex(Gbeta)$predictedSex)
table(pData(WB.noob)$Sex,getSex(Gbeta)$predictedSex)

#phenodata
pheno <- as.data.frame(cbind(Sex=pData(WB.noob)$Sex, Plate_ID=pData(WB.noob)$Plate_ID, cellprop))
# 1 female, 2 male

#quick check of the distribution of gender between plates
counts <- table(pheno[,"Sex"], pheno[,"Plate_ID"])
Percentage <- prop.table(counts, 2); 
barplot(Percentage, main="Distribution of sex within plates",
        xlab="plate", col=c("grey","white"), 
        legend= c("F","M"),args.legend = list(x = "topleft"));

pheno[,"Sex"]<-ifelse(pheno[,"Sex"]==2,0,1)
# 1 female, 0 male

## Cleaning up the methylation data
#' Filters a matrix of beta values by distance to SNP. Also removes crosshybridising probes and sex-chromosome probes.
betas.clean<-rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG <- dim(betas.clean)[1]
nCpG

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on sex (as a predictor of methylation)

#' first we can run a linear regression on a single cpg
j <- 133211

CpG.level <- betas.clean[j,]
CpG.name <- rownames(betas.clean)[j]
CpG.name

#' difference in methylation between different Sex
#' some statistics
knitr::kable(cbind(Min=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],min)),3),
      Mean=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],mean)),3), 
      Median=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],median)),3),
      Max=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],max)),3),
      SD=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],sd)),3),
      N=table(pheno[,"Sex"])))

boxplot(CpG.level ~ pheno[,"Sex"])
summary(lm(CpG.level~pheno[,"Sex"]))

#' EWAS and results 

#' sex as predictor  
#' note that CpGassoc is quite fast for running almost a half million regressions!
system.time(results1 <- cpg.assoc(betas.clean,pheno[,"Sex"]))

#' look at the results
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3]))
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3])[order(results1$results[,3]),])
#' check with previous results
cbind(results1$coefficients[j,4:5], results1$results[j,c(1,3)])

#' Bonferroni significant hits
table(results1$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results1$results[,5] < 0.05)


#' EWAS with adjustment for cell types
results2 <- cpg.assoc(betas.clean,pheno[,"Sex"], 
                    covariates=pheno[,"CD8T"]+ pheno[,"CD4T"] +  
                      pheno[,"NK"]  + pheno[,"Bcell"] + 
                      pheno[,"Mono"] + pheno[,"Gran"] + pheno[,"nRBC"])

#' look at the results
head(cbind(results2$coefficients[,4:5], P.value=results2$results[,3])[order(results2$results[,3]),])
#' Bonferroni significant hits
table(results2$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results2$results[,5] < 0.05)

#' qqplot and lambda interpretation  
#+ fig.width=13, fig.height=7, dpi=300
par(mfrow=c(1,1))
#plot(results1, main="QQ plot for association between methylation and sex")
plot(results2, main="QQ plot for association between methylation and sex \n adjusted for cells proportion")

#' Lambda - this is a summary of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline?  
#' estimated at -log10(0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)
#' Lambda for the first EWAS
lambda(results1$results[,3])
#' Lambda after the cell type adjustment
lambda(results2$results[,3])

# Volcano Plot-results2
plot(results2$coefficients[,4],-log10(results2$results[,3]), 
     xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot- adj for CellProp")
#Bonferroni treshold
abline(h = -log10(0.05/(nCpG)), lty=1, col="red", lwd=2)
#FDR treshold
abline(h = -log10(max(results1$results[results1$results[,5] < 0.05,3])), lty=1, col="blue", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
#' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
IlluminaAnnot = data.frame(
  chr=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations$chr,
  pos=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations$pos,
  Relation_to_Island=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Islands.UCSC$Relation_to_Island,
  UCSC_RefGene_Name=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$UCSC_RefGene_Name,
  UCSC_RefGene_Group=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$UCSC_RefGene_Group,
  probeA=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$ProbeSeqA,
  probeB=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$ProbeSeqB,
  Forward_Sequence=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$Forward_Sequence,
  SourceSeq=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$SourceSeq)

## Create CpG name and annotate row names
rownames(IlluminaAnnot) <- rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest)
IlluminaAnnot$Name <-rownames(IlluminaAnnot)
dim(IlluminaAnnot)

IlluminaAnnot = IlluminaAnnot [intersect(rownames(IlluminaAnnot), rownames(betas.clean)),]
dim(IlluminaAnnot)
datamanhat <- data.frame(CpG=results2$results[,1],Chr=as.character(IlluminaAnnot$chr),
                     Mapinfo=IlluminaAnnot$pos,Pval=results2$results[,3])
#' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))

#' Manhattan plot 
manhattan(datamanhat,"Chr","Mapinfo", "Pval", "CpG", main="Manhattan Plot - adj for CellProp")

#' End of script 02
#' 