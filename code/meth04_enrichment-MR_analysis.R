#' ---
#' title: Going Beyond EWAS with Enrichment Analyses and Causal Inference
#' author: Dr. Cavin Ward-Caviness
#' date: June 27, 2017
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' <style>
#' .myTable td {
#'  padding-left: 20px;
#'  padding-right: 20px;
#'  padding-top: 2px;
#'  padding-bottom: 2px;
#'  border:1px solid black;
#'  text-align: center;
#'  vertical-align: center;
#' }
#' 
#' .myTable tr:nth-child(even) {
#' background-color: #f2f2f2
#' } 
#' 
#' .myTable th {
#' background-color: black;
#'   color: white;
#'  padding-left: 20px;
#'  padding-right: 20px;
#'  padding-top: 2px;
#'  padding-bottom: 2px;
#'  border:1px solid white;
#'  text-align: center;
#' }  
#' </style>
  
  
#' This code performs.  
#' - Comparison of overlap of 450k and 850k features   
#' - Enrichment analyses for previous CHARGE EWAS  
#' - Enrichment analyses for smoking GWAS ([NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/))  
#' - Enrichment analyses for genomic feature  

# Install Packages for MR
if (!require(jsonlite)) install.packages("jsonlite")  # Installs devtools package if not already installed         
library(jsonlite)
if (!require(devtools)) install.packages("devtools")  # Installs devtools package if not already installed        
library(devtools)

install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
install_github("MRCIEU/MRInstruments")
library(MRInstruments)

#' ## Check for EWAS results and load necessary libraries
if(!exists("results2")){
  source("meth02_analyze_data.R")
}
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(shinyMethyl)
  library(ENmix)
  library(missMethyl) ### for pathway analyses
  library(CpGassoc) 
  library(lmtest)
  library(sva)
  library(limma)
  
  source("methylation-enrichment-EPIC.R") ### code adapted from Roby Joehanes for EPIC array
})

#' ## Read in necessary data
#'Published smoking EWAS data from [Epigenetic Signatures of Cigarette Smoking](https://doi.org/10.1161/CIRCGENETICS.116.001506 )

chrg_fmr_nvr <- read.table(file="../data/gwas_ewas_data/smk_former_vs_never.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, fill=TRUE)

#' GWAS data downloaded from the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) on 4/10/2017
smk_gwas_loc <- "../data/gwas_ewas_data/gwas-association-downloaded_2017-04-10-smoking.tsv"

#' Get 450K annotation
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))
IlluminaAnnot$Methyl450_Loci[IlluminaAnnot$Methyl450_Loci==TRUE] <- "450K"
IlluminaAnnot$Methyl450_Loci[IlluminaAnnot$Methyl450_Loci!="450K"] <- "EPIC only"

#' For ease I just directly merge annotation with the results data
results.anno <- merge(results2$results, IlluminaAnnot, by.x="CPG.Labels", by.y="Name")


#' ## Gene Ontology (GO) pathway enrichment.  
#' Here we will analyze the "suggestive" loci because there are only a few FDR loci making pathway enrichments suspect
results.sug <- subset(results.anno, P.value < 1E-5)
#' First this is done without controlling for the number of CpGs in each gene.  
#' Note: This produces a biased and somewhat inflated result.  
results.go.enrich.bias <- gometh(sig.cpg=results.sug$CPG.Labels, all.cpg=IlluminaAnnot$Name, collection="GO", array.type="EPIC", prior.prob=FALSE)
knitr::kable(head(topGO(results.go.enrich.bias)), format = "html", table.attr='class="myTable"')
#'    
#' Now control for the known bias of some genes having more annotated CpGs than others
results.go.enrich <- gometh(sig.cpg=results.sug$CPG.Labels, all.cpg=IlluminaAnnot$Name, collection="GO", array.type="EPIC", prior.prob=TRUE)
#+ results='asis'
knitr::kable(head(topGO(results.go.enrich)), format = "html", table.attr='class="myTable"')
#'    
#' Now to examine enrichment of different genomic features
#' First lets just examine the distribution of 450K vs EPIC only loci by CpGi annotation
# All CpGs

genome.locs <- table(results.anno$Methyl450_Loci, results.anno$Relation_to_Island)
knitr::kable(round(genome.locs/rowSums(genome.locs)*100,1), format = "html", table.attr='class="myTable"')
#' Table of all EPIC CpGs by genomic location and whether they were on the 450K array.  
#' Table cells represent percentages

#'
genome.locs <- table(results.sug$Methyl450_Loci, results.sug$Relation_to_Island)
#+ results='asis'
knitr::kable(round(genome.locs/rowSums(genome.locs)*100,1), format = "html", table.attr='class="myTable"')
#' Table of Suggestive (P < 1E-5) associations by genomic location and whether they were on the 450K array.  
#' Table cells represent percentages
#'
#' ## Genomic Feature Enrichment
#' ### Ehancer regions with two-sided test
#' Enhancer annotations are from the FANTOM 5 group
IlluminaAnnot$P5_YesNo <- ifelse(IlluminaAnnot$Phantom5_Enhancers=="",NA,TRUE)

# add column for CpGs with P < 1E-5
IlluminaAnnot$SMK_Sug <- ifelse(IlluminaAnnot$Name%in%results.sug$CPG.Labels, TRUE, FALSE)
# make 2x2 table
P5.enhancer_tbl <- with(IlluminaAnnot, table(SMK_Sug, !is.na(P5_YesNo), useNA = "ifany"))

#' Two sided p-value test for enrichment and depletion
doublemidp.test(P5.enhancer_tbl)

#' ### All regions one sided tests
#' Below comment is long running code which can be do beforehand and simply loaded
met_annot <- prep_annotation(IlluminaAnnot)

#' Code for one sided test is in additional R functions loaded at beginning of scrpit
cpgfeature_enrichment <- cpg_enrichment(results.sug$CPG.Labels, met_annot)
cpgfeature_enrichment <- cpgfeature_enrichment[order(cpgfeature_enrichment$Enrichment_P),]
#+ results='asis'
knitr::kable(cpgfeature_enrichment[c(1:10),], format = "html", table.attr='class="myTable"')

#' ## Compare with previous smoking EWAS
#' The data frame chrg_fmr_nvr only contains the FDR significant results from the former vs never smoking EWAS.  
#' Will now merge those with the associations from our EWAS with P < 0.05.  
results.fn.overlap <- merge(subset(results.anno, P.value < 0.05), chrg_fmr_nvr, by.x="CPG.Labels", by.y="Probe.ID")

#' ### How many results are FDR significant in both
sum(results.fn.overlap$FDR.x < 0.05) ### Number FDR significant in both
results.fn.overlap$CPG.Labels[results.fn.overlap$FDR.x < 0.05]

#' ### How many of the former associations do we replicate at a nominal (P < 0.05) level?
dim(results.fn.overlap)[1]

#' ### Results specific to this analysis
n.fn.overlap <- subset(results.sug, !CPG.Labels%in%c(chrg_fmr_nvr$Probe.ID))

#' Previous results were 450K only, so consider where in the genome those were as compared to current associations
overlap.genome <- table(results.fn.overlap$Methyl450_Loci, results.fn.overlap$Relation_to_Island)
n.overlap.genome <- table(n.fn.overlap$Methyl450_Loci, n.fn.overlap$Relation_to_Island)

#' Table of associations by genomic location that were FDR significant in CHARGE & have P < 0.05 in our smaller EWAS.  
#' Results given as percentages
#+ results='asis'
knitr::kable(round(overlap.genome/rowSums(overlap.genome)*100,1), format = "html", table.attr='class="myTable"')

#' "Suggestive" (P < 1E-5) associations from our EWAS that were not observed in the CHARGE EWAS by genomic location.  
#' Results given as percentages
#+ results='asis'
knitr::kable(round(n.overlap.genome/rowSums(n.overlap.genome)*100,1), format = "html", table.attr='class="myTable"')

#' ### Replicated associations from previous EWAS
#' Can check the number of replicated results.  
#' Use a bonferroni level of 0.05/# CHARGE FDR results
sum(results.fn.overlap$P.value.x < 0.05/dim(chrg_fmr_nvr)[1])

#' ## Overlap with smoking GWAS
gwas <- prep_gwas(read.delim(smk_gwas_loc, header=TRUE, sep="\t", as.is=TRUE));
gwas_enrichment <- cpg_gwas_enrichment(results.sug$CPG.Labels, met_annot, gwas)
#' ### Overall GWAS enrichment
#+ results='asis'
knitr::kable((attr(gwas_enrichment, "overall_enrichment_p")))
#' ### Top enrichment results
#' Was only one enrichment observed so only printing off the top row
#+ results='asis'
knitr::kable(gwas_enrichment[1,c("Diseases","N_Genes_In_GWAS","N_Genes_In_Dataset","Enrichment_P")], format = "html", table.attr='class="myTable"')


#' ## Mendelian Randomization (MR) Analysis
#' Conducted using the same code as used in [MR-base](http://www.mrbase.org/).  
#' We are using the R code version of MR-base found in the [TwoSampleMR Github](https://github.com/MRCIEU/TwoSampleMR).  
#' Hypothesis: Smoking has a causal effect on DNA methylation.  
#' MR Model: Smoking Variant (Instrument)  -> Smoking (Exposure) -> Methylation (Outcome).  


#'    
#'     
#' Load gwas catalog data to get IVs.  
#' Data can be downloaded and subsetted to the phenotype of interest, see commented code. This was pre-done here.  
# data(gwas_catalog)
# smk.gwas <- subset(gwas_catalog, grepl("smoking",Phenotype))
load("../data/MR_data/MRBase_GWAS-Catalog_Smoking.RData")

#' ### Step 1: Format the GWAS data as the exposure
init.gwas <- format_data(subset(smk.gwas, Phenotype=="Smoking behavior (smoking initiation)"))
cess.gwas <- format_data(subset(smk.gwas, Phenotype=="Smoking behavior (smoking cessation)"))

#' ### Step 2: Load the (pre-processed & re-analyzed) meQTL data as the outcome 
#' Data taken from the [BIOS eQTL browser](http://www.genenetwork.nl/biosqtlbrowser/)
# Below code will give an error. Is only because the effect allele frequency is not included. Can be ignored.
outcome.dat <- read_outcome_data(
  snps=c(init.gwas$SNP,cess.gwas$SNP),
  ## This file is a pre-formatted file of SNPs that are associated with methylation loci restricted to SNPs from smoking GWAS
  filename="../data/MR_data/Meta_Processed_BIOS_Init_Cess_GWAS_meQTLs.txt",
  sep="\t",
  snp_col="SNP",
  beta_col="Beta",
  se = "SE",
  effect_allele_col="effect_allele",
  pval_col="Pvalue",
  other_allele_col="other_allele"
)

#' ### Step 3: Restrict to CpGs with P < 0.05 in the EWAS, Harmonize the data, and perform the MR analysis
#' Need to show some level of association to be considered for a causal effect
outcome.dat <- subset(outcome.dat, outcome%in%results.anno$CPG.Labels[results.anno$P.value < 0.05])
dim(outcome.dat)[1]

#' Harmonize the data.  
#' Harmonization involves allele and strand matching & run the analysis on smoking cessation and initiation
if(outcome.dat$SNP%in%cess.gwas$SNP)
{
  
  cess.dat <- harmonise_data(
    exposure_dat=cess.gwas,
    outcome_dat = outcome.dat
  )
  # Run the MR analysis
  cess.MR <- mr(cess.dat)
  cess.MR <- merge(results.anno[,c("CPG.Labels","T.statistic","P.value","UCSC_RefGene_Name")], cess.MR, by.x="CPG.Labels", by.y="outcome")
}
#' With MR you must always check the output to make sure that effect direction from the original association (EWAS here) the direction of the causal effect.  
#' Smoking variants were taken from a GWAS was for smoking cessation while our EWAS was for non-smokers vs smokers so opposite effect directions (comparing columns T.statistic and b) are in fact consistent due to the "flipped" phenotypes.  
#' "T.statistic" is the effect of smoking on methylation estimated from our EWAS whereas "b"  is the causal effect estimated from the Mendelian Randomization. 
#' The units of T.statistic and b are NOT the same so the magnitude of effect cannot be compared, but the direction can and we do observe a causal effect of smoking on methylation.
#' The two rows represent tso similar methods to test for a causal effect.  
#' Information on MR and the difference in the methods is published in [Haycock PC, et al. (2016) Am J Clin Nutr. 103(4):965-978](http://ajcn.nutrition.org/content/103/4/965.short)
#+ results='asis'
knitr::kable(cess.MR, format = "html", table.attr='class="myTable"')
#' If we used multiple variants could generate a report describing several sensitivity analyses.  
#' `mr_report(cess.MR, output_path=".", output_type="html", author="Me", study="Smoking and methylation")`  

#' The only SNP retained was a cessation smoking variant so this will get skipped
# if(outcome.dat$SNP%in%init.dat$SNP)
# {
#  init.dat <- harmonise_data(
#    exposure_dat=init.gwas,
#    outcome_dat = outcome.dat
#  )
#  
#  init.MR <- mr(init.dat)
#  init.MR <- merge(results.anno[,c("CPG.Labels","T.statistic","P.value","UCSC_RefGene_Name")], init.MR, by.x="CPG.Labels", by.y="outcome")
#}

