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
model = model.matrix( ~smoker+sex+CD4+CD8+NK+B+MO+GR,data=pheno)


#' Introduction to limma 
#' see [Smyth GK. Stat Appl Genet Mol Biol 2004](https://www.ncbi.nlm.nih.gov/pubmed/16646809).  
suppressMessages(library(limma,minfi))
EWAS.limma <- eBayes(lmFit(betas.clean, design=model))
topTable(EWAS.limma, coef=2, number=Inf, sort.by="p")[1:10,]


#'Regions are now agglomerated from groups of significant probes 
#'Let's run the regional analysis using the Beta-values from our preprocessed data
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
#'Let's look at the second region
results.ranges[2]

# set up the grouping variables and colours
cols = c("magenta","red")[pheno$smoker]
names(cols) = levels(pheno$smoker)[pheno$smoker]

#'Draw the plot for a  DMR in\
#+ fig.width=8, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=2, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=1:nrow(pheno))

#'Draw the plot for another DMR\
#+ fig.width=8, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=3, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=1:nrow(pheno))

#' cleanup
rm(tx.hg19,tx.hg38,tx.mm10,snpsall,myBetas,myannotation,crosshyb,XY.probes);gc()


#'Extracting CpGs-names and locations
coord = dmrcoutput.smoking$results$coord[1]
coord = stri_match(coord,regex="^(chr.+):(\\d+)-(\\d+)$")

chr = coord[2]
start = as.integer(coord[3])
end = as.integer(coord[4])

#'CpG ID and individual metrics

cpgs = subset(dmrcoutput.smoking$input, CHR == chr & pos >= start & pos <= end)
knitr::kable(cpgs)


#' ## Predicting smoking with EpiSmokEr: https://www.biorxiv.org/content/10.1101/487975v1.article-info
# Make sure rows of pheno match betas column names
rownames(pheno)<-pheno$gsm
identical(colnames(beta),rownames(pheno))

# pheno needs a column for sex,in the format of 1 and 2 representing men and women respectively
pheno$sex<-ifelse(pheno$sex=="m",1,2)
# 121 CpGs are used selected by LASSO along with Sex to get 3 categories (current, former and never smokers)
result <- epismoker(dataset=beta, samplesheet = pheno, method = "SSt")
result[,]
# Let's look how well the prediction performed
table(pheno$smoker,result$PredictedSmokingStatus)


#' End of script 03