.libPaths("C:/EBC4/Rpackages")

#+ setdir01, echo = F
knitr::opts_knit$set(root.dir = "../")

options(warn=-1)
library(data.table)
suppressMessages(library(limma))
suppressMessages(library(ENmix))
options(warn=0)

#' Load the data
load("data/processed.rda")

#'# Exploring global DNA Methylation variability via PCs
# Let's look at the effect of sex
#+ fig.width=8, fig.height=6, dpi=300
plotMDS(beta, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pheno$sex)],
        dim=c(1,2),cex=1.5)
legend("center", legend=levels(factor(pheno$sex)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))

# Let's look at the effect of sex but let's remove the X and Y chromosomes
betas.clean = beta[manifest[probe_type=="cg" & !chr %in% c("X","Y")]$index,]
plotMDS(betas.clean, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pheno$sex)],
        dim=c(1,2),cex=1.5)
legend("center", legend=levels(factor(pheno$sex)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))

#' It would be useful to look at several traits with global variability
cov<-data.frame(pheno[,2:9])
npc <- 20 # Top 20 PCs
svd <- prcomp(t(na.omit(betas.clean)))
screeplot(svd, npc, type = "barplot")
eigenvalue <- svd[["sdev"]]^2
prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue)) * 100
cat("Top ", npc, " principal components can explain ", 
    prop, "% of data \n    variation", "\n")
screeplot(svd,npc,type="barplot")
pcrplot(na.omit(betas.clean), cov, npc=5)


#'# Epigenetic Age
#'Load package for Age-Prediction
library(wateRmelon)
DNAmAge<-as.vector(agep(beta))
hist(DNAmAge)
boxplot(DNAmAge);stripchart(DNAmAge, vertical = T,method = "jitter", add = T, pch = 20, col = 'red')

#'Agreement between Horvath's Epigenetic age and Hannum's clock
data(hannumCoef)
length(hannumCoef)
DNAmAge.Hannum<-as.vector(agep(beta,coeff=hannumCoef,method = "hannum"))


#' Correlation; agreement
cor.test(DNAmAge,DNAmAge.Hannum)
plot(DNAmAge.Hannum,DNAmAge,pch=21,ylab="Hortvath's DNAm Age",
     xlab="Hannum's DNAm Age",cex=1.2, bg=alpha("deepskyblue",0.45),main="Epigenetic Clocks")
legend("topleft",legend=paste0("r=",round(cor(DNAmAge.Hannum,DNAmAge),2)),bty="n")
abline(lm(DNAmAge~DNAmAge.Hannum),col="red",lw=2)


#' Age Acceleration Residuals
AgeAccelerationResidual<-residuals(lm(DNAmAge.Hannum~DNAmAge))
boxplot(AgeAccelerationResidual ~pheno$smoker, col=c("blue","red"))
wilcox.test(AgeAccelerationResidual ~ pheno$smoker)

#' Differences by Smoking status
boxplot(DNAmAge ~pheno$smoker, col=c("blue","red"))
wilcox.test(DNAmAge ~ pheno$smoker)
boxplot(DNAmAge.Hannum ~pheno$smoker, col=c("blue","red"))
wilcox.test(DNAmAge.Hannum ~ pheno$smoker)
