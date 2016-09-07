#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  
#+ setdir02, echo = F
knitr::opts_knit$set(root.dir = "../")

#' we have a processed dataset with 15 samples (otherwise we run script 01)
if(!exists("WB.noob")){
  source("code/meth01_process_data.R")
}
dim(WB.noob)

#' load packages
suppressPackageStartupMessages({
  library(CpGassoc)
  library(data.table)
  library(qqman)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(DMRcate)
  library(MASS) 
  library(sandwich) 
  library(lmtest) 
})

#'## predict sex from methylation
Gbeta <- mapToGenome(WB.noob)  #map to the genome
#' getSex predicts sex based on X and Y chromosome methylation intensity
getSex(Gbeta) 
#' we see that our predictions match the phenodata
table(pData(WB.noob)$Sex,getSex(Gbeta)$predictedSex)
# cleanup
rm(Gbeta)

#' consolidate our phenodata
pheno <- as.data.frame(cbind(Sex=pData(WB.noob)$Sex, Plate_ID=pData(WB.noob)$Plate_ID, cellprop))
pheno[,"Sex"] <- ifelse(as.numeric(pheno[,"Sex"])==2,0,1)
# 1 female, 0 male

# cleanup
rm(WB.noob)

#' quick check of the distribution of gender between plates
table(pheno[,"Sex"], pheno[,"Plate_ID"])

## Cleaning up the methylation data
#' Filters a matrix of beta values by distance to SNP. Also removes crosshybridising probes and sex-chromosome probes.
dim(betas.rcp)
betas.clean <- rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG <- dim(betas.clean)[1]
nCpG

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on sex (as a predictor of methylation)  

#' First we can run a linear regression on a single CpG that we have already picked
j <- 133211
CpG.level <- betas.clean[j,]
CpG.name <- rownames(betas.clean)[j]
CpG.name

#' difference in methylation between males and females for this CpG
#' some descriptive statistics
knitr::kable(cbind(Min=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],min)),3),
                   Mean=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],mean)),3), 
                   Median=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],median)),3),
                   Max=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],max)),3),
                   SD=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],sd)),3),
                   N=table(pheno[,"Sex"])))
#' difference in beta methylation values between different Sex  
boxplot(CpG.level ~ pheno[,"Sex"], main="Beta-values")

#' linear regression on betas
summary(lm(CpG.level~pheno[,"Sex"]))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' comparison with m-values
CpG.mlevel <- log2(betas.clean[j,])-log2(1-betas.clean[j,])

knitr::kable(cbind(Min=round(simplify2array(tapply(CpG.mlevel, pheno[,"Sex"],min)),3),
                   Mean=round(simplify2array(tapply(CpG.mlevel, pheno[,"Sex"],mean)),3), 
                   Median=round(simplify2array(tapply(CpG.mlevel, pheno[,"Sex"],median)),3),
                   Max=round(simplify2array(tapply(CpG.mlevel, pheno[,"Sex"],max)),3),
                   SD=round(simplify2array(tapply(CpG.mlevel, pheno[,"Sex"],sd)),3),
                   N=table(pheno[,"Sex"])))

boxplot(CpG.mlevel ~ pheno[,"Sex"], main="M-values")

#' linear regression on m-values
summary(lm(CpG.mlevel~pheno[,"Sex"]))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' we can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
summary(lm(CpG.level~pheno[,"Sex"]))$adj.r.squared

#' model on mvalues
summary(lm(CpG.mlevel~pheno[,"Sex"]))$adj.r.squared

#'## EWAS and results using CpGassoc
#'see [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  

#' sex as predictor  
#' note that CpGassoc is quite fast for running almost a half million regressions!
system.time(results1 <- cpg.assoc(betas.clean, pheno[,"Sex"]))

#' there are several components of the results
class(results1)
names(results1)
#' look at a few results  
#' here effect size is ~ mean difference in methylation proportion
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3]))
#' and the top hits
head(cbind(results1$coefficients[,4:5], P.value=results1$results[,3])[order(results1$results[,3]),])
#' check with previous result on our selected CpG (running lm without CpGassoc)
cbind(results1$coefficients[j,4:5], results1$results[j,c(1,3)])
summary(lm(CpG.level~pheno[,"Sex"]))

#' Bonferroni significant hits
table(results1$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results1$results[,5] < 0.05)

#' EWAS with adjustment for cell types
#' now we can run the linear regression on betas adjusting for cell proportions
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
#' we can also see them with:
results2$FDR.sig

#' results
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
IlluminaAnnot = data.frame(
  chr=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations$chr,
  pos=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations$pos,
  Relation_to_Island=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Islands.UCSC$Relation_to_Island,
  UCSC_RefGene_Name=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$UCSC_RefGene_Name,
  UCSC_RefGene_Group=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other$UCSC_RefGene_Group)

## Create CpG name and annotate row names
rownames(IlluminaAnnot) <- rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest)
IlluminaAnnot$Name <-rownames(IlluminaAnnot)
dim(IlluminaAnnot)

IlluminaAnnot = IlluminaAnnot [intersect(rownames(IlluminaAnnot), rownames(betas.clean)),]
dim(IlluminaAnnot)
datamanhat <- data.frame(CpG=results2$results[,1],Chr=as.character(IlluminaAnnot$chr),
                         Mapinfo=IlluminaAnnot$pos, UCSC_RefGene_Name=IlluminaAnnot$UCSC_RefGene_Name, 
                         Pval=results2$results[,3], Eff.Size = results2$coefficients[,4], Std.Error = results2$coefficients[,5])

#' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))

#' see where the top hits are
head(datamanhat[order(datamanhat$Pval), ])

#' ## Genomic inflation in EWAS
#' qqplot and lambda interpretation  
#+ fig.width=13, fig.height=7, dpi=300
par(mfrow=c(1,1))
plot(results1, main="QQ plot for association between methylation and sex")
plot(results2, main="QQ plot for association between methylation and sex \n adjusted for cell proportions")

#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline?  
#' estimated at -log10(0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)
#' Lambda for the first EWAS
lambda(results1$results[,3])
#' Lambda after cell type adjustment
lambda(results2$results[,3])

#' Volcano Plot-results2
#' with Bonferroni threshold and current FDR
plot(results2$coefficients[,4],-log10(results2$results[,3]), 
     xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot \n adjusted for cell proportions")
#Bonferroni threshold & FDR threshold
abline(h = -log10(0.05/(nCpG)), lty=1, col="#FDE725FF", lwd=2)
abline(h = -log10(max(results2$results[results2$results[,5] < 0.05,3])), lty=1, col="#440154FF", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
#' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
manhattan(datamanhat,"Chr","Mapinfo", "Pval", "CpG", 
          suggestiveline = -log10(max(results2$results[results2$results[,5] < 0.05,3])), 
          genomewideline = -log10(0.05/(nCpG)), 
          main = "Manhattan Plot \n adjusted for cell proportions")

#' cleanup
rm(j, nCpG, CpG.name, CpG.level, CpG.rlm, CpG.mlevel, lm.fit.rob.bayes, datamanhat, IlluminaAnnot, IlluminaHumanMethylation450kanno.ilmn12.hg19, lambda)
#' End of script 02
#'  