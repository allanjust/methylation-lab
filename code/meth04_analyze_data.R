#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  
#+ setdir02, echo = F
#' Local library
.libPaths("C:/EBC4/Rpackages")
knitr::opts_knit$set(root.dir = "../")

#' we have a processed dataset with 30 samples (otherwise we run script 01)
#if(!exists("data/processed.rds")){
#  source("code/meth01_process_data.R")
#}

# load the data
load("data/processed.rda")

#' load packages
suppressPackageStartupMessages({
  library(DMRcate) # for regional analysis
  library(magrittr)
  library(CpGassoc) # for running association analysis between methylation levels values and phenotype of interest
  library(data.table) # for fast aggregation of large data 
  library(qqman) # for visualization of data
  library(stringi) # string manipulation
})

#' Code categorical variable as factors
pheno$smoker %<>% factor
pheno$sex    %<>% factor

## Cleaning up the methylation data
#' Filters a matrix of beta values by distance to single nucleotide polymorphism (SNP) and SNPs with the minor allele frequency (MAF) of 5% (rare variant). 
#' Also removes crosshybridising probes and sex-chromosome probes.
dim(beta)
betas.clean = beta[manifest[probe_type=="cg" & !chr %in% c("X","Y")]$index,]
nCpG = dim(betas.clean)[1]
nCpG

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on Smoking status (as a predictor of methylation)  

#' First we can run a linear regression on a single CpG that we have already picked
CpG.name = "cg05575921"
pheno$CpG.level <- betas.clean[CpG.name,]

#' Difference in methylation between smokers and non-smokers for this CpG
#' some descriptive statistics

pheno[,.(
   Min    = min   (CpG.level) %>% round(3)
  ,Mean   = mean  (CpG.level) %>% round(3)
  ,Median = median(CpG.level) %>% round(3)
  ,Max    = max   (CpG.level) %>% round(3)
  ,SD     = sd    (CpG.level) %>% round(3)
  ,.N),by=smoker] %>% knitr::kable(.)

#' Difference in beta methylation values between Smokers and non smokers
boxplot(CpG.level ~ smoker,pheno,main=paste0("Beta-values\n", CpG.name), col=c("blue","red"),ylim=c(.3,1))

#' Linear regression on betas
lm(CpG.level ~ smoker,data=pheno) %>% summary %>% coef

#' Comparison with m-values
pheno[,CpG.mlevel:=log2(CpG.level/(1-CpG.level))]

pheno[,.(
   Min    = min   (CpG.mlevel) %>% round(3)
  ,Mean   = mean  (CpG.mlevel) %>% round(3)
  ,Median = median(CpG.mlevel) %>% round(3)
  ,Max    = max   (CpG.mlevel) %>% round(3)
  ,SD     = sd    (CpG.mlevel) %>% round(3)
  ,.N),by=smoker] %>% knitr::kable(.)


par(mfrow=c(1,2))
boxplot(CpG.level  ~ smoker,data=pheno,main=paste0("Beta-values\n",CpG.name), col=c("blue","red"))
boxplot(CpG.mlevel ~ smoker,data=pheno,main=paste0("M-values\n"   ,CpG.name), col=c("blue","red"))

#' linear regression on m-values
lm(CpG.mlevel ~ smoker,data=pheno) %>% summary %>% coef

#' We can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
lm(CpG.level  ~ smoker,data=pheno) %>% summary %$% adj.r.squared

#' model on mvalues
lm(CpG.mlevel ~ smoker,data=pheno) %>% summary %$% adj.r.squared

#'## EWAS and results using CpGassoc
#'see [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  

#' Smoking as predictor  
#' note that CpGassoc is quite fast for running almost a million regressions!

pheno[,smoke_dummy:=ifelse(smoker=="smoker",1,0)]

system.time(results1 <- cpg.assoc(betas.clean, pheno$smoke_dummy,fdr.cutoff=0.1))

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
summary(lm(CpG.level~smoker,pheno))

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
          ,pheno$smoke_dummy
          ,covariates=as.data.frame(pheno[,.(sex,CD8,CD4,NK,B,MO,GR)])
          )

print(results2)

#'using mvalues
results3 <- cpg.assoc(
           betas.clean
          ,pheno$smoke_dummy
          ,covariates=as.data.frame(pheno[,.(sex,CD8,CD4,NK,B,MO,GR)])
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


datamanhat = cbind(results2$results,results2$coefficients)
setDT(datamanhat)
datamanhat = datamanhat[,.(probe_id=CPG.Labels,effect.size,std.error,P.value)]

datamanhat = merge(datamanhat,manifest[,.(probe_id,chr,mapinfo)],by="probe_id")

#' See where the top hits are
datamanhat[order(P.value)][1:7]

#' Volcano Plot-results2
#' Bonferroni threshold
#' 
plot(-log10(P.value) ~ effect.size,data=datamanhat,xlab="Estimate",ylab="-log10(p-value)",main="Volcano Plot\nadjusted for cell proportions",ylim=c(0,8))
abline(h = -log10(0.05/(nCpG)), lty=1, col="#FDE725FF", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
#' Cast the variable chr (so we can simplify and use a numeric x-axis)
datamanhat[,chr:=as.integer(chr)]

qqman::manhattan(datamanhat,chr="chr",bp="mapinfo",p="P.value",snp="probe_id"
   ,suggestiveline=FALSE, genomewideline = -log10(0.05/(nCpG)),ylim=c(0,8)
   ,main = "Manhattan Plot \n adjusted for cell proportions")

#' cleanup
rm(nCpG,CpG.name,datamanhat,lambda,results1,results2,results3)
gc()
#' End of script 02
