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
# To Fix lirary IlluminaHumanMethylation450kanno.ilmn12.hg19
suppressPackageStartupMessages({
  library(CpGassoc)
  library(data.table)
  library(qqman)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(glmnet)
  library(DMRcate)
  library(MASS) 
  library(sandwich) 
  library(lmtest) 
})


#' consolidate our phenodata
pheno <- as.data.frame(cbind(Smoke=pData(WB.noob)$SMOKE_STATUS, Sex=pData(WB.noob)$SEX, Age=pData(WB.noob)$AGE,
                             Array=pData(WB.noob)$Array, cellprop))
pheno$Smoke <- ifelse(pheno$Smoke=="SMOKER",1,0)
# 1 Smoker , 0 non-smokers

# cleanup
rm(WB.noob)

#' quick check of the distribution of gender between plates
table(pheno[,"Smoke"], pheno[,"Array"])

## Cleaning up the methylation data
#' Filters a matrix of beta values by distance to SNP/MAF. Also removes crosshybridising probes and sex-chromosome probes.
dim(betas.rcp)
betas.clean <- rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG <- dim(betas.clean)[1]
nCpG

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on Smoking status (as a predictor of methylation)  

#' First we can run a linear regression on a single CpG that we have already picked
j <- 19094
CpG.level <- betas.clean[j,]
CpG.name <- rownames(betas.clean)[j]
CpG.name

#' difference in methylation between males and females for this CpG
#' some descriptive statistics
knitr::kable(cbind(Min=round(simplify2array(tapply(CpG.level, pheno[,"Smoke"],min)),3),
                   Mean=round(simplify2array(tapply(CpG.level, pheno[,"Smoke"],mean)),3), 
                   Median=round(simplify2array(tapply(CpG.level, pheno[,"Smoke"],median)),3),
                   Max=round(simplify2array(tapply(CpG.level, pheno[,"Smoke"],max)),3),
                   SD=round(simplify2array(tapply(CpG.level, pheno[,"Smoke"],sd)),3),
                   N=table(pheno[,"Smoke"])))
#' difference in beta methylation values between Smokers and non smokers
boxplot(CpG.level ~ pheno[,"Smoke"], main="Beta-values", col=c("blue","red"), xaxt="n")
axis(1,at=c(1,2),adj=1,labels=cbind("Non-smoker","Smoker"))

#' linear regression on betas
summary(lm(CpG.level~pheno[,"Smoke"]))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' comparison with m-values
CpG.mlevel <- log2(betas.clean[j,])-log2(1-betas.clean[j,])

knitr::kable(cbind(Min=round(simplify2array(tapply(CpG.mlevel, pheno[,"Smoke"],min)),3),
                   Mean=round(simplify2array(tapply(CpG.mlevel, pheno[,"Smoke"],mean)),3), 
                   Median=round(simplify2array(tapply(CpG.mlevel, pheno[,"Smoke"],median)),3),
                   Max=round(simplify2array(tapply(CpG.mlevel, pheno[,"Smoke"],max)),3),
                   SD=round(simplify2array(tapply(CpG.mlevel, pheno[,"Smoke"],sd)),3),
                   N=table(pheno[,"Smoke"])))

boxplot(CpG.mlevel ~ pheno[,"Smoke"], main="M-values")

#' linear regression on m-values
summary(lm(CpG.mlevel~pheno[,"Smoke"]))$coefficients[2,c("Estimate", "Pr(>|t|)","Std. Error")]

#' we can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
summary(lm(CpG.level~pheno[,"Smoke"]))$adj.r.squared

#' model on mvalues
summary(lm(CpG.mlevel~pheno[,"Smoke"]))$adj.r.squared

#'## EWAS and results using CpGassoc
#'see [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  

#' Smoking as predictor  
#' note that CpGassoc is quite fast for running almost a half million regressions!
system.time(results1 <- cpg.assoc(betas.clean, pheno[,"Smoke"]))

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
summary(lm(CpG.level~pheno[,"Smoke"]))

#' Bonferroni significant hits
table(results1$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results1$results[,5] < 0.05)

#' EWAS with adjustment for cell types
#' now we can run the linear regression on betas adjusting for cell proportions
#' Need sex as indicator in covariate matrix
pheno$Sex <- ifelse(pheno$Sex=="male",1,0)
results2 <- cpg.assoc(betas.clean,pheno[,"Smoke"], 
                        covariates=pheno[,"Sex"] + pheno[,"Age"]+
                        pheno[,"CD8T"]+ pheno[,"CD4T"] +
                        pheno[,"NK"]  + pheno[,"Bcell"] + 
                        pheno[,"Mono"] + pheno[,"Gran"])

#' look at the results
head(cbind(results2$coefficients[,4:5], P.value=results2$results[,3])[order(results2$results[,3]),])
#' Bonferroni significant hits
table(results2$results[,3] < 0.05/(nCpG))
#' FDR significant hits
table(results2$results[,5] < 0.05)
#' we can also see them with:
results2$FDR.sig

#'using mvalues
results3 <- cpg.assoc(betas.clean, pheno[,"Smoke"], 
                      covariates=pheno[,"Sex"] + pheno[,"Age"]+
                        pheno[,"CD8T"]+ pheno[,"CD4T"] +
                        pheno[,"NK"]  + pheno[,"Bcell"] + 
                        pheno[,"Mono"] + pheno[,"Gran"], 
                      logit.transform = TRUE)

#' look at the results
head(cbind(results3$coefficients[,4:5], P.value=results3$results[,3])[order(results3$results[,3]),])
#' Bonferroni significant hits
table(results3$results[,3] < 0.05/(nCpG))
#' FDR significant hits
results3$FDR.sig

#' results
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))

#' Restrict to good quality probes and order data frames
IlluminaAnnot <- IlluminaAnnot[IlluminaAnnot$Name %in% results2$results$CPG.Labels,]
IlluminaAnnot <- IlluminaAnnot[match(results2$results$CPG.Labels, IlluminaAnnot$Name),]

#' Check that CpGs are align
identical(IlluminaAnnot$Name,results2$results$CPG.Labels)

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
plot(results1, main="QQ plot for association between methylation and Smoking")
plot(results2, main="QQ plot for association between methylation and Smoking \n adjusted for cell proportions")
plot(results3, main="QQ plot for association between (mvals) methylation and Smoking \n adjusted for cell proportions")

#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline?  
#' estimated at -log10(0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)
#' Lambda for the first EWAS
lambda(results1$results[,3])
#' Lambda after cell type adjustment
lambda(results2$results[,3])
#' Lambda after cell type adjustment with mvalues
lambda(results3$results[,3])

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
rm(j, nCpG, CpG.name, CpG.level, CpG.rlm, CpG.mlevel, lm.fit.rob.bayes, datamanhat, IlluminaAnnot, IlluminaHumanMethylation450kanno.ilmn12.hg19, 
   lambda)
#' End of script 02
#'  
