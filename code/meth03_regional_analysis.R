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
#'  see [Peters et al. Bioinformatics 2015](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-6).  
#' Other popular options for conducting Regional DNA methylation analysis in R are Aclust and bumphunter 
suppressMessages(library(DMRcate)) # Popular package for regional DNA methylation analysis


#' First we need to define a model
model <- model.matrix(~as.factor(pheno$Smoke)+
                      as.factor(pheno$Sex)+
                      as.numeric(pheno$Age)+
                      as.numeric(pheno$CD8T)+
                      as.numeric(pheno$NK)+
                      as.numeric(pheno$Bcell)+
                      as.numeric(pheno$Mono)+
                      as.numeric(pheno$Gran))

#'Regions are now agglomerated from groups of significant probes 
#'Let's run the regional analysis using the Beta-values from our preprocessed data
myannotation <- cpg.annotate("array", betas.clean, analysis.type="differential",arraytype="EPIC",
                             what="Beta",design=model, coef=2)


#'We don't find any significant regions (FDR<0.05), So let's try a simpler model as an example
model <- model.matrix(~as.factor(pheno$Smoke)+
                      as.factor(pheno$Sex)+
                      as.numeric(pheno$Age))

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
pheno$Smoke <- ifelse(pheno$Smoke==1, "Smoker", "Non-Smoker")
# set up the grouping variables and colours
colors<-c("magenta","red")
groups <- colors[1:length(unique(pheno$Smoke))]
names(groups) <- levels(factor(pheno$Smoke))
cols <- groups[as.character(factor(pheno$Smoke))]
samps <- 1:nrow(pheno)
#'Draw the plot for the top DMR
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=samps)


#'Extracting CpGs-names and locations
chr <- gsub(":.*", "", dmrcoutput.smoking$results$coord[1])
start <- gsub("-.*", "", gsub(".*:", "", dmrcoutput.smoking$results$coord[1]))
end <- gsub(".*-", "", dmrcoutput.smoking$results$coord[1])
#'CpG ID and individual metrics
cpgs <- dmrcoutput.smoking$input[dmrcoutput.smoking$input$CHR %in% chr & dmrcoutput.smoking$input$pos >= start & dmrcoutput.smoking$input$pos <=end,]
knitr::kable(cpgs[1:4,])

#' Load package for regional analysis "Bumphunter"
#'  see [Jaffe et al. Int J Epidemiol. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22422453). 
suppressMessages(library("bumphunter","minfi","registerDoSEQ","doParallel"))
registerDoSEQ()
#'Create ratioset from clean betas
data.rs <- RatioSet(Beta = betas.clean,annotation=c(array= "IlluminaHumanMethylationEPIC",  annotation = "ilm10b2.hg19")); # create RatioSet                                                                                      
data.grs <- mapToGenome(data.rs); # create GenomicRatioSet  


#'Bumphunter using 14% DNAm difference
Cores<-detectCores()
registerDoParallel(cores = Cores-1)
dmrs.10 <- bumphunter(data.grs, design = model, coef=2,nullMethod="bootstrap",cutoff = 0.14, B=10, type="Beta") # 29 bumps at 20% difference in methylation
#'Look at top DMRs
dmrs.10$table[1:4,]


#' End of script 03
#' 