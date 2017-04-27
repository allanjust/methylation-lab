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
chrg_curr_nvr <- read.table(file="M:/Epigenetics Boot Camp/BootCampDemo/smk_curr_vs_never.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, fill=TRUE)
chrg_fmr_nvr <- read.table(file="M:/Epigenetics Boot Camp/BootCampDemo/smk_former_vs_never.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, fill=TRUE)
## GWAS data downloaded from the NHGRI-EBI GWAS Catalog on 4/10/2017
smk_gwas_loc <- "M:/Epigenetics Boot CAmp/BootCampDemo/gwas-association-downloaded_2017-04-10-smoking.tsv"

### assumes that the model is named "unadj.out" and is the ouput of a model from CpGassoc
### MethylSetEx.noob is the MethylSet for the data created after noob background correction
### noob.rcp is the matrix of betas used in the analysis, has the CpGs that were retained for the analysis

### get annotation
annot <- getAnnotation(MethylSetEx.noob)
annot$Methyl450_Loci[annot$Methyl450_Loci==TRUE] <- "450K"
annot$Methyl450_Loci[annot$Methyl450_Loci!="450K"] <- "EPIC only"

### merge with annotation and label which are 450k loci and which are EPIC specific
unadj.out <- merge(unadj.out, annot, by.x="CPG.Labels", by.y="Name")

### limit to just the "suggetive" loci - here P < 0.0005 just as an example
unadj.sug <- subset(unadj.out, P.value < 0.0005)


### perform GO pathway enrichment while controlling for the number of CpGs per gene
unadj.gst <- gometh(sig.cpg=unadj.sug$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC", plot.bias=TRUE)
topGO(unadj.gst)
### GO pathway enrichment without controlling for the number of CpGs per gene (biased result)
unadj.gst.bias <- gometh(sig.cpg=unadj.sug$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC", prior.prob=FALSE)
topGO(unadj.gst.bias)

### make frequency table of the genomic locations by which array the loci came from
genome.locs <- table(unadj.out$Methyl450_Loci, unadj.out$Relation_to_Island)
round(genome.locs/rowSums(genome.locs)*100,1)

### enhancer enrichment
annot$X450k_Enhancer[annot$X450k_Enhancer==""]<-NA
annot$P5_YesNo <- ifelse(annot$Phantom5_Enhancers=="",NA,TRUE)

with(annot, table(P5_YesNo, useNA = "ifany"))
with(annot, table(X450k_Enhancer, useNA = "ifany"))

with(annot, print(prop.table(table(
  P5_YesNo, useNA = "ifany")), digits = 2))

### add column for hits
annot$SMK_Sug <- ifelse(annot$Name%in%unadj.sug$CPG.Labels, TRUE, FALSE)
### make 2x2 table
P5.enhancer_tbl <- with(annot, table(SMK_Sug, !is.na(P5_YesNo), useNA = "ifany"))
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
unadj.fn.overlap <- merge(subset(unadj.out, P.value < 0.05), chrg_fmr_nvr, by.x="CPG.Labels", by.y="Probe.ID")
sum(unadj.fn.overlap$P.value < 0.05/dim(chrg_fmr_nvr)[1]) 

### can also see which associations are specific to this analysis
n.fn.overlap <- subset(unadj.sug, !CPG.Labels%in%c(chrg_fmr_nvr$Probe.ID))

### now remake frequency tables to see where in the genome the replicated and analysis-specific loci are
overlap.genome <- table(unadj.fn.overlap$Methyl450_Loci, unadj.fn.overlap$Relation_to_Island)
n.overlap.genome <- table(n.fn.overlap$Methyl450_Loci, n.fn.overlap$Relation_to_Island)
round(overlap.genome/rowSums(overlap.genome)*100,1)
round(n.overlap.genome/rowSums(n.overlap.genome)*100,1)

### GO enrichment of replicated sites from former vs never smoker EWAS
unadj.gst.fn <- gometh(sig.cpg=unadj.fn.overlap$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC")
topGO(unadj.gst.fn)

## now can do the same with the current vs never smoking GWAS from CHARGE which may better fit our phenotype
### which CHARGE results (curr vs never) do we replicate at a nominal level?
unadj.cn.overlap <- merge(subset(unadj.out, P.value < 0.05), chrg_curr_nvr, by.x="CPG.Labels", by.y="Probe.ID")
sum(unadj.cn.overlap$P.value < 0.05/dim(chrg_curr_nvr)[1]) ### none bonferroni

### which results are not in the current vs never smoking EWAS
n.cn.overlap <- subset(unadj.sug, !CPG.Labels%in%c(chrg_curr_nvr$Probe.ID))

### compare percent loci replicated of the two EWAS
rrate <- c(dim(unadj.fn.overlap)[1]/dim(chrg_fmr_nvr)[1], dim(unadj.cn.overlap)[1]/dim(chrg_curr_nvr)[1])*100
names(rrate) <- c("F v N", "C v N")

### genomic location frequency of overlapping loci for CHARGE current vs never smoking EWAS
overlap.genome <- table(unadj.cn.overlap$Methyl450_Loci, unadj.cn.overlap$Relation_to_Island)
n.cn.overlap.genome <- table(n.cn.overlap$Methyl450_Loci, n.cn.overlap$Relation_to_Island)
round(overlap.genome/rowSums(overlap.genome)*100,1)
round(n.cn.overlap.genome/rowSums(n.cn.overlap.genome)*100,1)

### GO pathway enrichment of replicated CpGs from EWAS of current vs never smokers
unadj.gst.cn <- gometh(sig.cpg=unadj.cn.overlap$CPG.Labels, all.cpg=rownames(noob.rcp), collection="GO", array.type="EPIC")
topGO(unadj.gst.cn)

#### enrichment for genes found in smoking GWAS (NHGRI)
gwas <- prep_gwas(read.delim(smk_gwas_loc, header=TRUE, sep="\t", as.is=TRUE));
gwas_enrichment <- cpg_gwas_enrichment(unadj.sug$CPG.Labels, met_annot, gwas)
print(attr(gwas_enrichment, "overall_enrichment_p"))



