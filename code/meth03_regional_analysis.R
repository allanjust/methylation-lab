#'# Regional DNA methylation analysis using DMRcate and bumphunter  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R & meth02_process_data.R 
#+ setdir03, echo = F
knitr::opts_knit$set(root.dir = "../")

#'  we have already set up our analysis
if(!exists("pheno")){
  source("code/meth02_analyze_data.R")
}

#' Load package for regional analysis "DMRcate"
#' see [Peters et al. Bioinformatics 2015](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-6).  
#' Other popular options for conducting Regional DNA methylation analysis in R are Aclust and bumphunter 
suppressMessages(library(DMRcate)) # Popular package for regional DNA methylation analysis


#' First we need to define a model
model = model.matrix( ~SMOKE_STATUS+SEX+AGE+CD8T+NK+Bcell+Mono+Gran,data=pheno )


#'Regions are now agglomerated from groups of significant probes 
#'Let's run the regional analysis using the Beta-values from our preprocessed data
myannotation <- cpg.annotate("array", betas.clean, analysis.type="differential",arraytype="EPIC",
                             what="Beta",design=model, coef=2)

#' Introduction to limma 
#' see [Smyth GK. Stat Appl Genet Mol Biol 2004](https://www.ncbi.nlm.nih.gov/pubmed/16646809).  
suppressMessages(library(limma,minfi))
EWAS.limma <- eBayes(lmFit(betas.clean, design=model))
topTable(EWAS.limma, coef=2, number=Inf, sort.by="p")[1:10,]

#'We don't find any significant regions (FDR<0.05), So let's try a simpler model as an example
model = model.matrix( ~SMOKE_STATUS+SEX+AGE,data=pheno )

EWAS.limma <- eBayes(lmFit(betas.clean, design=model))
topTable(EWAS.limma, coef=2, number=Inf, sort.by="p")[1:10,]

myannotation <- cpg.annotate("array", betas.clean, analysis.type="differential",arraytype="EPIC",
                             what="Beta",design=model, coef=2)

#'Regions are now agglomerated from groups of significant probes 
#'where the distance to the next consecutive probe is less than lambda nucleotides away
dmrcoutput.smoking <- dmrcate(myannotation, lambda=1000, C=2)

#'Let's look at the results
head(dmrcoutput.smoking$results)

#'Visualizing the data can help us understand where the region lies 
#'relative to promoters, CpGs islands or enhancers

#' Let's extract the genomic ranges and annotate to the genome
results.ranges <- extractRanges(dmrcoutput.smoking, genome = "hg19")

#' Plot the DMR using the Gviz

#' if you are interested in plotting genomic data the Gviz is extremely useful
#'Let's look at the first region
results.ranges[1]

# set up the grouping variables and colours
cols = c("magenta","red")[pheno$SMOKE_STATUS]
names(cols) = levels(pheno$SMOKE_STATUS)[pheno$SMOKE_STATUS]

#'Draw the plot for the top DMR
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=1:nrow(pheno))


#'Extracting CpGs-names and locations
coord = dmrcoutput.smoking$results$coord[1]
coord = stri_match(coord,regex="^(chr.+):(\\d+)-(\\d+)$")

chr = coord[2]
start = as.integer(coord[3])
end = as.integer(coord[4])

#'CpG ID and individual metrics

cpgs = subset(dmrcoutput.smoking$input, CHR == chr & pos >= start & pos <= end)
knitr::kable(cpgs)

#' Load package for regional analysis "Bumphunter"
#'  see [Jaffe et al. Int J Epidemiol. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22422453). 
suppressMessages(library("bumphunter","minfi","registerDoSEQ","doParallel"))
registerDoSEQ()
#'Create ratioset from clean betas
data.rs <- RatioSet(Beta = betas.clean,annotation=c(array= "IlluminaHumanMethylationEPIC",  annotation = "ilm10b2.hg19")); # create RatioSet                                                                                      
data.grs <- mapToGenome(data.rs); # create GenomicRatioSet  


#'Bumphunter using 14% DNAm difference
Cores<-detectCores()
dmrs.10 <- bumphunter(data.grs, design = model, coef=2,nullMethod="bootstrap",cutoff = 0.14, B=10, type="Beta") # 29 bumps at 20% difference in methylation
#'Look at top DMRs
dmrs.10$table[1:4,]

#'Load package for Age-Prediction
library(wateRmelon)
DNAmAge<-as.vector(agep(betas.rcp))

#'Agreement between chronological age and DNAm-Age
pheno$Age<-as.numeric(as.character(pheno$Age))
cbind(DNAmAge,pheno$Age,pheno$Smoke)

#' Correlation
cor.test(DNAmAge,as.numeric(as.character(pheno$Age)))
plot(pheno$Age,DNAmAge,pch=21,ylab="Age Predicted",
     xlab="Age Reported",cex=1.2, bg=alpha("deepskyblue",0.45),main="Prediction of Age")
legend("topleft",legend=c("r=",paste(round(cor(DNAmAge,pheno$Age),2))),bty="n")
abline(lm(DNAmAge~pheno$Age),col="red",lw=2)


#' Age Acceleration Residuals
AgeAccelerationResidual<-residuals(lm(DNAmAge~pheno$Age))
boxplot(AgeAccelerationResidual ~pheno$Smoke, col=c("red","blue"))
wilcox.test(AgeAccelerationResidual ~ pheno$Smoke)
t.test(AgeAccelerationResidual ~ pheno$Smoke)


#' Online Calculator
Horvath<-read.csv("C:/EBC3/methylation-lab/data/EpigeneticAge/MethylationData.output.csv")
cor.test(Horvath$DNAmAge,as.numeric(as.character(pheno$Age)))
plot(pheno$Age,Horvath$DNAmAge,pch=21,ylab="Age Predicted",
     xlab="Age Reported",cex=1.2, bg=alpha("deepskyblue",0.45),main="Prediction of Age")
legend("topleft",legend=c("r=",paste(round(cor(Horvath$DNAmAge,pheno$Age),2))),bty="n")
abline(lm(Horvath$DNAmAge~pheno$Age),col="red",lw=2)



#' End of script 03
#' 