## Comparison of overlap of 450k and 850k features
## enrichment analyses for previous CHARGE EWAS
## enrichment analyses for smoking GWAS (NHGRI-EBI GWAS Catalog 4/10/2017)
## enrichment analyses for genomic feature

suppressPackageStartupMessages({
  library(EPICdemo)
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(shinyMethyl)
  library(ENmix)
  library(FlowSorted.Blood.450k)
  library(missMethyl) ### for pathway analyses
  library(CpGassoc) 
  library(lmtest)
  library(sva)
  library(limma)
  library(cate)
  library(data.table)
  
  source("M:/Epigenetics Boot Camp/BootCampDemo/methylation-enrichment-EPIC.R") ### code adapted from Roby Joehanes for EPIC array
})

### read in necessary data
## smoking data from "Epigenetic Signatures of Cigarette Smoking" (https://doi.org/10.1161/CIRCGENETICS.116.001506 )
chrg_fmr_nvr <- read.table(file="M:/Epigenetics Boot Camp/BootCampDemo/smk_former_vs_never.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, fill=TRUE)
## GWAS data downloaded from the NHGRI-EBI GWAS Catalog on 4/10/2017
smk_gwas_loc <- "M:/Epigenetics Boot CAmp/BootCampDemo/gwas-association-downloaded_2017-04-10-smoking.tsv"

### assumes that the model is named "results2" and is the ouput of a model from CpGassoc
### Gbeta is the GenomicRatioSet for the data created after noob background correction

if(!exists("results2")){
  source("code/meth02_analyze_data.R")
}

### get annotation
IlluminaAnnot<-as.data.frame(getAnnotation(Gbeta))
IlluminaAnnot$Methyl450_Loci[IlluminaAnnot$Methyl450_Loci==TRUE] <- "450K"
IlluminaAnnot$Methyl450_Loci[IlluminaAnnot$Methyl450_Loci!="450K"] <- "EPIC only"

### merge with annotation and label which are 450k loci and which are EPIC specific
results.anno <- merge(results2$results, IlluminaAnnot, by.x="CPG.Labels", by.y="Name")
results.anno$Methyl450_Loci[results.anno$Methyl450_Loci==TRUE] <- "450K"
results.anno$Methyl450_Loci[results.anno$Methyl450_Loci!="450K"] <- "EPIC only"

### limit to just the "suggetive" loci - here P < 1E-5 just as an example
results.sug <- subset(results.anno, P.value < 1E-5)


### perform GO pathway enrichment while controlling for the number of CpGs per gene
results.go.enrich <- gometh(sig.cpg=results.sug$CPG.Labels, all.cpg=results.anno$CPG.Labels, 
                            collection="GO", array.type="EPIC", plot.bias=TRUE)
topGO(results.go.enrich)
### GO pathway enrichment without controlling for the number of CpGs per gene (biased result)
results.go.enrich.bias <- gometh(sig.cpg=unadj.sug$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC", prior.prob=FALSE)
topGO(results.go.enrich.bias)

### make frequency table of the genomic locations by which array the loci came from
genome.locs <- table(results.anno$Methyl450_Loci, results.anno$Relation_to_Island)
round(genome.locs/rowSums(genome.locs)*100,1)

### enhancer enrichment
IlluminaAnnot$X450k_Enhancer[IlluminaAnnot$X450k_Enhancer==""]<-NA
IlluminaAnnot$P5_YesNo <- ifelse(IlluminaAnnot$Phantom5_Enhancers=="",NA,TRUE)

with(IlluminaAnnot, table(P5_YesNo, useNA = "ifany"))
with(IlluminaAnnot, table(X450k_Enhancer, useNA = "ifany"))

with(annot, print(prop.table(table(
  P5_YesNo, useNA = "ifany")), digits = 2))

### add column for hits
IlluminaAnnot$SMK_Sug <- ifelse(IlluminaAnnot$Name%in%results.sug$CPG.Labels, TRUE, FALSE)
### make 2x2 table
P5.enhancer_tbl <- with(IlluminaAnnot, table(SMK_Sug, !is.na(P5_YesNo), useNA = "ifany"))
addmargins(P5.enhancer_tbl)
### print proportion of hits
print(prop.table(P5.enhancer_tbl, 1), digits = 2)

### two sided p-value test for enrichment and depletion
doublemidp.test(P5.enhancer_tbl)

### one sided test for enrichment and depletion (separately) enrichment of suggestive CpGs for genomic features
met_annot <- prep_annotation(annot)
cpgfeature_enrichment <- cpg_enrichment(unadj.sug$CPG.Labels, met_annot)

### compare to CHARGE former vs never smoking EWAS
### which results do we replicate at a nominal (P < 0.05) level?
results.fn.overlap <- merge(subset(results.anno, P.value < 0.05), chrg_fmr_nvr, by.x="CPG.Labels", by.y="Probe.ID")
dim(results.fn.overlap)[1]
sum(results.fn.overlap$FDR.x < 0.05) ### Number FDR significant in both

### can also see which associations are specific to this analysis
n.fn.overlap <- subset(results.sug, !CPG.Labels%in%c(chrg_fmr_nvr$Probe.ID))

### now remake frequency tables to see where in the genome the replicated and analysis-specific loci are
overlap.genome <- table(results.fn.overlap$Methyl450_Loci, results.fn.overlap$Relation_to_Island)
n.overlap.genome <- table(n.fn.overlap$Methyl450_Loci, n.fn.overlap$Relation_to_Island)
round(overlap.genome/rowSums(overlap.genome)*100,1)
round(n.overlap.genome/rowSums(n.overlap.genome)*100,1)

### GO enrichment of replicated sites from former vs never smoker EWAS
results.go.enrich.fn <- gometh(sig.cpg=results.fn.overlap$CPG.Labels, all.cpg=results.anno$CPG.Labels, collection="GO", array.type="EPIC")
topGO(results.go.enrich.fn)

rrate <- dim(results.fn.overlap)[1]/dim(chrg_fmr_nvr)[1]*100
paste0(round(rrate,1),"%")

#### enrichment for genes found in smoking GWAS (NHGRI)
gwas <- prep_gwas(read.delim(smk_gwas_loc, header=TRUE, sep="\t", as.is=TRUE));
gwas_enrichment <- cpg_gwas_enrichment(results.sug$CPG.Labels, met_annot, gwas)
print(attr(gwas_enrichment, "overall_enrichment_p"))
gwas_enrichment[c(1:5),c("Diseases","N_Genes_In_GWAS","N_Genes_In_Dataset","Enrichment_P")]

## now can do the same with the current vs never smoking GWAS from CHARGE which may better fit our phenotype
### which CHARGE results (curr vs never) do we replicate at a nominal level?
# chrg_curr_nvr <- read.table(file="M:/Epigenetics Boot Camp/BootCampDemo/smk_curr_vs_never.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, fill=TRUE)
#unadj.cn.overlap <- merge(subset(unadj.out, P.value < 0.05), chrg_curr_nvr, by.x="CPG.Labels", by.y="Probe.ID")
#sum(unadj.cn.overlap$P.value < 0.05/dim(chrg_curr_nvr)[1]) ### none bonferroni
#
### which results are not in the current vs never smoking EWAS
#n.cn.overlap <- subset(unadj.sug, !CPG.Labels%in%c(chrg_curr_nvr$Probe.ID))
#
### compare percent loci replicated of the two EWAS
#rrate <- c(dim(results.fn.overlap)[1]/dim(chrg_fmr_nvr)[1], dim(unadj.cn.overlap)[1]/dim(chrg_curr_nvr)[1])*100
#names(rrate) <- c("F v N", "C v N")
#
### genomic location frequency of overlapping loci for CHARGE current vs never smoking EWAS
#overlap.genome <- table(unadj.cn.overlap$Methyl450_Loci, unadj.cn.overlap$Relation_to_Island)
#n.cn.overlap.genome <- table(n.cn.overlap$Methyl450_Loci, n.cn.overlap$Relation_to_Island)
#round(overlap.genome/rowSums(overlap.genome)*100,1)
#round(n.cn.overlap.genome/rowSums(n.cn.overlap.genome)*100,1)

### GO pathway enrichment of replicated CpGs from EWAS of current vs never smokers
#results.go.enrich.cn <- gometh(sig.cpg=unadj.cn.overlap$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC")
#topGO(results.go.enrich.cn)

###########################################
#
# MR analyses using MR-base
# there is an online version (www.mrbase.org)
# but we will be using the R code available via GitHub
#
# Hypothesis: Smoking has a causal effect on DNA methylation
# MR Model: IV  -> Exposure -> Outcome
#           SNP -> Smoking  -> Methylation
##########################################

# step 1: Install necessary libraries and code from GitHub
if (!require(jsonlite)) install.packages("jsonlite")          # Installs devtools package if not already installed
library(jsonlite)
if (!require(devtools)) install.packages("devtools")          # Installs devtools package if not already installed
library(devtools)

install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
install_github("MRCIEU/MRInstruments")
library(MRInstruments)

# step 2: download gwas catalog data to get IVs
data(gwas_catalog)
smk.gwas <- subset(gwas_catalog, grepl("smoking",Phenotype))
## alternative
## load(data/MRBase_GWAS-Catalog_Smoking.RData)

# step 3: format the GWAS data as the exposure
init.gwas <- format_data(subset(smk.gwas, Phenotype=="Smoking behavior (smoking initiation)"))
cess.gwas <- format_data(subset(smk.gwas, Phenotype=="Smoking behavior (smoking cessation)"))

# step 4: load the (pre-processed & re-analyzed) meQTL data as the outcome 
## data taken from the BIOS eQTL browser: http://www.genenetwork.nl/biosqtlbrowser/
outcome.dat <- read_outcome_data(
  snps=c(init.gwas$SNP,cess.gwas$SNP),
  filename="M:/Epigenetics Boot Camp/BootCampDemo/Meta_Processed_BIOS_Init_Cess_GWAS_meQTLs.txt",
  sep="\t",
  snp_col="SNP",
  beta_col="Beta",
  se = "SE",
  effect_allele_col="effect_allele",
  pval_col="Pvalue",
  other_allele_col="other_allele"
)

## step 5: restrict to CpGs with a nominal (P < 0.05) signal in the EWAS
nom_only=TRUE
if(nom_only)
{
  outcome.dat <- subset(outcome.dat, outcome%in%results.anno$CPG.Labels[results.anno$P.value < 0.05])
}

### step 6: hamonise the data (make sure alleles match) & run the analysis on smoking cessation and initiation
if(outcome.dat$SNP%in%cess.gwas$SNP)
{
  cess.dat <- harmonise_data(
    exposure_dat=cess.gwas,
    outcome_dat = outcome.dat
  )
  cess.MR <- mr(cess.dat)
  cess.MR <- merge(results.anno[,c("CPG.Labels","T.statistic","P.value","UCSC_RefGene_Name")], cess.MR, by.x="CPG.Labels", by.y="outcome")
  
  ### if we used multiple variants could generate a report describing several sensitivity analyses
  ## mr_report(cess.MR, output_path=".", output_type="html", author="Me", study="Smoking and methylation")
}

if(outcome.dat$SNP%in%init.dat$SNP)
{
  init.dat <- harmonise_data(
    exposure_dat=init.gwas,
    outcome_dat = outcome.dat
  )
  
  init.MR <- mr(init.dat)
  init.MR <- merge(results.anno[,c("CPG.Labels","T.statistic","P.value","UCSC_RefGene_Name")], init.MR, by.x="CPG.Labels", by.y="outcome")
}

## check output to make sure that effect direction from EWAS matches that from the MR
cess.MR


