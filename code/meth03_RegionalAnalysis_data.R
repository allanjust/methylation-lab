#'# Regional DNA methylation analysis  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  & meth02_process_data.R 
#'  
rmarkdown::render("meth02_analyze_data.R")

#' Load package for regional analysis "DMRcate"
#'  see Peters et al. Bioinformatics 2015. https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-6
#' Other popular options for conducting Regional DNA methylation analysis are Aclust and bumphunter 
suppressMessages(library(DMRcate)) # Popular package for regional DNA methylation analysis


#' First we need to define a model
model<-model.matrix(~as.factor(pheno$Sex)+
                     as.numeric(pheno$CD8T)+
                     as.numeric(pheno$NK)+
                     as.numeric(pheno$Bcell)+
                     as.numeric(pheno$Mono)+
                     as.numeric(pheno$Gran)+
                     as.numeric(pheno$nRBC))

#'Let's run the regional analysis
myannotation <- cpg.annotate("array", Mvals.ComBat, analysis.type="differential",
                             design=model, coef=2)

#'Regions are now agglomerated from groups of significant probes where the distance to the next consecutive probe is less than lambda nucleotides away
dmrcoutput.sex <- dmrcate(myannotation, lambda=1000, C=2)

#'Let's look at our results
head(dmrcoutput.sex$results)

#'Visualizing the data can help us understand where the region lies relative to promoters, CpGs islands or enhancers

#' Let's extract the genomic ranges and annotate to the genome
results.ranges <- extractRanges(dmrcoutput.sex, genome = "hg19")
results.ranges[2]

#' Plot the DMR using the Gviz

#'If you are interested in plotting genomic data the Gviz is extremely useful
#'Let's look at the first region
results.ranges[1]
pheno$Sex<-ifelse(pheno$Sex==1,"Female","Male")
groups <- c(Female="magenta", Male="forestgreen")
cols <- groups[as.character(pheno$Sex)]
#+ fig.width=9, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas.clean, phen.col=cols, genome="hg19")

#'Let's look at the second region
results.ranges[2]
#+ fig.width=9, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=2, CpGs=betas.clean, phen.col=cols, genome="hg19")
#' End of script
#' 