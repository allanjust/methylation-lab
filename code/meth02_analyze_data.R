#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  

#Elena Colicino

#' recall that we have a processed dataset with 15 samples
dim(WB)

suppressPackageStartupMessages({
  library(CpGassoc)
  library(data.table)
})

# # Evaluate any remaining batch effects from design using the principal component analysis
# PCobject = prcomp(na.omit(t(betas.bmiq)), retx = T, center = T, scale. = T);
# # Extract the Principal Components from SVD
# PCs = PCobject$x;
# 
# # Extract the proportion of variability explained by the top R=2 PC.
# R = 2
# propvar = summary(PCobject)$importance["Proportion of Variance", 1:R]
# sum(propvar)*100
# cummvar = summary(PCobject)$importance["Cumulative Proportion", 1:R]
# cummvar
# 
# X = 1  # First PC
# Y = 2  # Second PC
# 
# # Save Graphic (Plate Effects)
# variable = as.factor(WB$Plate_ID)
# myColors = rainbow(length(table(variable)))
# names(myColors) = levels(variable)
# # X-Y
# par(mfrow=c(1,1))
# plot(PCs[,X], PCs[,Y], xlab = paste("PC", X,"(",round(propvar[X]*100,1),"%)", sep =""),
#      ylab = paste("PC", Y, "(",round(propvar[Y]*100,1), "%)", sep =""), pch=16,
#      col = myColors)
# legend("topleft", c("Plate 1","Plate 2"), pch=c(21,21),bty='n', cex=1,pt.bg = myColors, col=myColors)

## Extract sex prediction from methylation data
# use median total intensity of the X chromosome-mapped probes and
# the median total intensity of the Y-chromosome-mapped probes

Gbeta<-mapToGenome(WB)  #map to the genome
estSex <- getSex(Gbeta) #get sex
estSex
Gbeta <- addSex(Gbeta, sex = estSex)  #add predicted sex to the phenodata
#note: predicted sex is already included in the dataset-> warning for replacment of that column

Gbeta$predictedSex<-as.factor(Gbeta$predictedSex)
cbind(PSex=Gbeta$predictedSex,Sex=Gbeta$Sex)

#phenodata
pheno<-cbind(Sex=Gbeta$predictedSex, Plate_ID=Gbeta$Plate_ID, cellprop)
# 1 female, 2 male

#quick check of the distribution of sex among plates
#check if sex is well distributed between plates
counts <- table(pheno[,"Sex"], pheno[,"Plate_ID"])
Percentage <- prop.table(counts, 2); 
barplot(Percentage, main="Distribution of sex within plates",
        xlab="plate", col=c("grey","white"), 
        legend= c("F","M"),args.legend = list(x = "topleft"));

pheno[,"Sex"]<-ifelse(pheno[,"Sex"]==2,0,1)
# 1 female, 0 male

#'# Running an Epigenome Wide Association
#' Here we run an EWAS on sex (as a predictor of methylation)

# first we can run a single regression on a single cpg
j<-2835

CpG.level<-betas.bmiq[j,]
CpG.name<-rownames(betas.bmiq)[j]

#difference in methylation between different Sex
rbind(Min=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],min)),3),
      Mean=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],mean)),3), 
      Median=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],median)),3),
      Max=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],max)),3),
      SD=round(simplify2array(tapply(CpG.level, pheno[,"Sex"],sd)),3),
      N=table(pheno[,"Sex"]))

boxplot(CpG.level ~ pheno[,"Sex"])
# plot(pheno[,"Sex"],CpG.level, main="Methylation and Sex", ylab="DNA methytilation level", xlab="Sex") 
# abline(lm(CpG.level~pheno[,"Sex"]), col="red")
summary(lm(CpG.level~pheno[,"Sex"]))

#EWAS 
results1<-cpg.assoc(betas.bmiq,pheno[,"Sex"])
results1

#check with previous results
cbind(results1$coefficients[j,4:5], results$results[j,c(1,3)])

#EWAS with adjustment for cell types
results2<-cpg.assoc(betas.bmiq,pheno[,"Sex"], 
                    covariates=pheno[,"CD8T"]+ pheno[,"CD4T"] +  pheno[,"NK"]  + pheno[,"Bcell"] + pheno[,"Mono"] + pheno[,"Gran"] + pheno[,"nRBC"])
results2

#' looking at results
 
Coef<-results1$coefficients[,4:5]
setnames(Coef, colnames(Coef), c("Est", "SD"))
Coef<-Coef[order(rownames(Coef)),]
Pval<-results$results[,c(1,3)]
setnames(Pval, colnames(Pval), c("CpG", "Pval"))
Pval<-Pval[order(Pval$CpG),]
identical(rownames(Coef), Pval$CpG)
#TRUE

Results_Sex<-data.frame(cbind(Coef,Pval))
save(Results_Sex, file=paste0(path.to.output,"Results_Sex_", Sys.Date(),".RData"))
write.csv(Results_Sex, file=paste0(path.to.output,"Results_Sex_", Sys.Date(),".csv"))



#' qqplot and lambda interpretation

# Lambda
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)
lambda(results$results[,3])

# Volcano Plot
nCpG=dim(B)[1]
plot(Results_Sex$Est,-log10(Results_Sex$Pval), ylim=c(0,5), xlab="Estimate", ylab="-log10(Pvalue)", main="Volcano Plot")
#canonical thresold
abline(h=-log10(0.05/(nCpG)), lty=1, col="red", lwd=2)
#relaxed thresold (just for us)
abline(h=-log10(0.05/(nCpG)), lty=1, col="blue", lwd=2)


#' volcano plot interpretation

#' looking at a top hit



#' End of script