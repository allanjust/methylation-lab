#' ---
#' title: "Processing, analyzing, and interpreting epigenome-wide DNA methylation data"
#' author: "Andrea Baccarelli, Andres Cardena, Elena Colicino, Jonathan Heiss, Allan Just"
#' date: "June, 2019"
#' geometry: margin=2cm
#' number_sections: true
#' ---
#' Local library
.libPaths("C:/EBC4/Rpackages")
#+ setdir01, echo = F
knitr::opts_knit$set(root.dir = "../")


options(warn=-1)
library(stringi)
library(magrittr)
library(data.table)
library(svd)
library(ewastools)
options(warn=0)


#' ## Importing the data
#' 1. Read in the file `data/pheno.csv` using `fread` from the data.table package, save it as object named `pheno`.
#' 2. Import the methylation data using the function `read_idats`.
#'    You only have to provide the first part of the file name without the `_Red.idat.gz` or `_Grn.idat.gz` suffixes.
#' 3.  Save as object named `meth`

pheno = fread("data/pheno.csv")
meth = read_idats("data/" %s+% pheno$gsm,quiet=TRUE)


#' Take a look at the dataset

# The name of the platform (450K/EPIC)
meth$platform

#' A manifest with probe IDs, color channel, genomic coordinates and other important information
meth$manifest[4000:4010]
#' Please note that Type I probes have `addressU` and `addressM` (two different beads), whereas Type II probes only have a single `addressU`

#' Most of the probes are of Type II
table(meth$manifest$channel)

#' Not all probes are targeting CpG sites
table(meth$manifest$probe_type)

#' Similar manifest for the control probes
head(meth$controls)

#' Also included
# Matrices contained fluorescence intensities for the methylated (`M`) and unmethylated (`U`) signals
dim(meth$M)
meth$U[201:203,1:3]
meth$M[201:203,1:3]

#' Matrices with the bead copy number (`N` and `V`)
dim(meth$N)
meth$N[201:203,1:3]

#' For some beads the copy number is zero. Because of the random assembly some did not end up on the chips
colSums(meth$N==0)

#' Matrices with the out-of-band intensities
dim(meth$oobG$M) # same number as Type I red probes
dim(meth$oobR$U)

#' Some meta data
head(meth$meta)

#' # Quality control

#' ## Control metrics
#' The first quality check evaluates 17 control metrics which are describe in the [BeadArray Controls Reporter Software Guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf) from Illumina. Exemplary, the "Bisulfite Conversion II" metric is plotted below. Three samples fall below the Illumina-recommended cut-off of 1. Input for `control_metrics()` is the output of `read_idats()`, e.g. the object holding raw or dye-bias-corrected intensities.

#' 1. Use the functions `control_metrics` and `sample_failure`
#' 2. Are there any failed samples?

meth %>% control_metrics %>% sample_failure -> pheno$failed
# Alternative syntax:
# pheno$failed = sample_failure(control_metrics(meth))

table(pheno$failed,useNA='ifany')
#' There are no failed samples.
#'

#' ## Sex check
#'
#' 1. Apply the function `check_sex` on the `meth` object to compute the normalized average total fluorescence intensity of probes targeting the X and Y chromosomes.
#' 2. Use the function `predict_sex` to infer the sex of the sample donors from the methylation data.
#' 3. Add column `predicted_sex` to `pheno`
#' 4. Are there samples where `sex != predicted_sex`? What is the sample ID?
#' 5. Flag the problematic samples with a logical variable `exclude` in `pheno`

pheno[,c("X","Y"):=check_sex(meth)]

pheno$predicted_sex = predict_sex(pheno$X,pheno$Y,which(pheno$sex=="m"),which(pheno$sex=="f"))

pheno[sex!=predicted_sex]
plot(Y~X,data=pheno,type="n")
text(Y~X,labels=sex,col=ifelse(sex=="m",2,1),data=pheno)

pheno[,exclude:=FALSE]

pheno[sex!=predicted_sex,.(gsm,sex,predicted_sex)]
pheno[sex!=predicted_sex,exclude:=TRUE] # flag sample

#' ## Detection p-values
#' 

meth = ewastools::detectionP(meth)
chrY = meth$manifest[chr=='Y',index]
detP = meth$detP[chrY,]
detP = colSums(detP<0.01,na.rm=TRUE)

boxplot(split(detP,pheno$predicted_sex),ylab="# of detected Y chromosome probes")
split(detP,pheno$predicted_sex) %>% sapply(mean)

#' Almost all of the 416 chromosome probes are called detected in male samples, for female samples on average 60 are called detected.

#' How many probes are undetected (not counting the Y chromosome). The cut-off used here is 0.01
meth$detP[-chrY,] %>% is_weakly_greater_than(0.01) %>% table(useNA="ifany")

#' Less than 0.2% are undetected
round((33749/(33749+16939303))*100,3)

#' We should mask these undetected probes.
meth = ewastools::mask(meth,0.01)

#' ## Dye-bias correction
#' 
#' Infinium BeadChips use two fluorescent dyes that are linked to the nucleotides used in the the single-base extension step. A and T nucleotides use are linked with a red dye (the red color channel), G and C nucleotides are linked with a green dye (green color channel). Uncorrected data usually feature higher intensities in the red color channel, the so-called dye bias. For probes of Infinium type II design, which use separate color channels to measure the methylated and unmethylated signal, this results in a shifted distribution of beta-values. (Probes of Infinium design type I are not affected, as they measure both signals in the same color channel.) Dye-bias correction normalizes the red and green color channel. `ewastools` provides an improved version of RELIC ([Xu et al., 2017](https://doi.org/10.1186/s12864-016-3426-3)) using robust Theil-Sen estimators.


color_bias = meth %>% dont_normalize
beta       = meth %>% correct_dye_bias %>% dont_normalize

#' We can look at a few methylation values on the fly and see whether dye-bias correction changed them
meth$manifest$channel[201:203] # One probe for each type/color channel
color_bias[201:203,1:3] %>% round(4)
beta      [201:203,1:3] %>% round(4)

#' If we calculate beta-values from raw data, we can observe the dye bias as a deviation of the beta-values for heterozygous SNPs from 0.5
snps = meth$manifest[probe_type=="rs" & channel=="Both"]$index

plot (density(color_bias[snps,14],na.rm=TRUE,bw=0.1),col=1,main="Dye-bias correction")
lines(density(beta      [snps,14],na.rm=TRUE,bw=0.1),col=2)
abline(v=0.5,lty=3)
legend("topleft",col=1:2,legend=c("raw","corrected"),lwd=1)

#' 
plot (density(beta[meth$manifest$channel=="Grn" ,1],na.rm=TRUE),col="green",main="Distribution of beta-values")
lines(density(beta[meth$manifest$channel=="Red" ,1],na.rm=TRUE),col="red")
lines(density(beta[meth$manifest$channel=="Both",1],na.rm=TRUE),col="black")
legend("topright",legend=c("Type I Red","Type I Grn","Type II"),lwd=1,col=c("red","green","black"))

#' ## SNP outliers
#' `snp_outliers()` returns the average log odds of belonging to the outlier component across all SNP probes. I recommend to flag samples with a score greater than -4 for exclusion.
snps = meth$manifest[probe_type=="rs"]$index
genotypes = call_genotypes(beta[snps,],learn=FALSE)
pheno$outlier = snp_outliers(genotypes)

stripchart(pheno$outlier,method="jitter",pch=4)
abline(v=-4,lty="dotted",col=2)

pheno[outlier>-4]
pheno[outlier>-4,exclude:=TRUE]

#' The one sample failing this check is the same sample that did not belong to either the male or female cluster in the plot above. This is strong evidence that this sample is indeed contaminated. While SNP outliers can also result from poorly performing assays, the sample passed the first quality check looking at the control metrics, therefore rendering this possibility unlikely. Another cause for a high outlier score is sample degradation (e.g., FFPE samples).

pheno$donor_id = enumerate_sample_donors(genotypes)
pheno[,n:=.N,by=donor_id]
pheno[n>1,.(gsm,donor_id)]

pheno[gsm=="GSM2260543",exclude:=TRUE] # drop duplicate

#' ## Principal component analysis
#' Principal component analysis is a popular feature reduction method: it projects high-dimensional data into a lower-dimensional representation while trying to retain as much variability as possible. This is especially useful when either individual features are highly correlated and it is therefore reasonable to summarize them, or when (sometimes subtle) traces of background effects can be found across of large number of features.
#' 1. Get of subset of beta without probes on the X and Y chromosome

chrXY = meth$manifest$chr %in% c("X","Y") | meth$manifest$probe_type == "rs"
pcs = beta[-chrXY,]
pcs = pcs - rowMeans(pcs)
pcs = na.omit(pcs)
pcs = t(pcs)
pcs = trlan.svd(pcs,neig=2)
pcs = pcs$u

pheno$pc1 = pcs[,1]
pheno$pc2 = pcs[,2]

plot(pc2~pc1,col=ifelse(sex=="m",2,1),data=pheno)
text(pc2~pc1,labels=pheno[34]$gsm,data=pheno[34],pos=4,offset=1,col=2)

pheno[gsm=="GSM2219539",exclude:=TRUE]

#' GSM2219539 is actually a lung tissue sample from another GEO dataset. It dominates the first principal component and should be excluded as it otherwise could drastically change the results of downstream analyses.
#' PCA may be applied iteratively. After excluding samples that manifest as outliers, repeating PCA can give very different principal components.
#'

#' ## Leukocyte composition
#' This quality check will only apply in case of blood samples (blood is, however, one of the most commonly studied tissues). The function `estimateLC()` implements the [Houseman method](https://doi.org/10.1186/1471-2105-13-86) to predict the leukocyte composition. The user has the choice between various sets of model parameters trained on various reference datasets (see `?estimateLC` for a list of options). The function operates on the matrix of beta-values.
#' 
#' 1. Estimate the leukocyte composition using the function `estimateLC`. What does the function return?
#' 2. Add the cell proportion estimates as new columns to the `pheno` data.table using the function `cbind`
#' 3. Plot the proportions of granulocytes versus (the column named `GR`) using
#' 4. Set the `exclude` flag TRUE for the conspicuous sample

LC = estimateLC(beta,ref="Reinius")

stripchart(rowSums(LC),xlab="Sum of cell proportions",m="jitter")
abline(v=1,lty=3,col=2)

pheno = cbind(pheno,LC)

plot(pheno$GR,ylim=c(0,1))

pheno[which.min(GR),.(gsm,exclude)]
#' This is the lung tissue sample from before

pheno[which.max(GR),.(gsm,exclude)]
#' This is actually a sample of purified granulocytes

pheno[gsm=="GSM1185585",exclude:=TRUE]


#' LC stratified by smoking status
#' 
LC$smoker = pheno$smoker
LC = melt(LC,value.name="proportion",variable.name="cell_type",id.vars="smoker")

boxplot(proportion ~ smoker+cell_type,LC,col=1:2,main="Cell type distribution by smoking status",xaxt="n")
axis(1,at=seq(from=1.5, to=11.5,by=2),adj=1,labels=unique(LC$cell_type))
legend("topright",c("Non-smoker","Smoker"),pch=15,bty='n',col=1:2)

# drop problematic samples
keep = which(!pheno$exclude)

# drop columns no longer needed
pheno = pheno[,.(gsm,smoker,sex,CD4,CD8,NK,MO,GR,B)]

# write data to disk
pheno = pheno[keep]
beta  = beta[,keep]
manifest = copy(meth$manifest)

#save(pheno,manifest,beta,file="data/processed.rda")
