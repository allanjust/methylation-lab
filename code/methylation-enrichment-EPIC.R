#############################################################
# Enrichment analysis for methylation
# Copyright (C) 2013-2017 Roby Joehanes
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, a copy is available at
# https://www.gnu.org/licenses/gpl-3.0.en.html
#
# Code modified by Cavin Ward-Caviness to work with EPIC data
#
#############################################################

# Assumed input files:
# met-annot.txt is Illumina methylation annotation file
# gwas_catalog_* is NHGRI/EBI GWAS catalog downloaded from the web

# Note: This is specifically for 450K chip. Others may need adjustments.

# Typical usage:

#library(data.table);
#source("methylation-enrichment-EPIC.R");
#met_annot <- prep_annotation(data.frame(fread("met-annot.txt")));
#gwas <- prep_gwas(read.delim("gwas_catalog_v1.0-downloaded_2015-11-03.tsv", header=TRUE, sep="\t", as.is=TRUE));
#
#tt <- cpg_gwas_enrichment(cpgs, met_annot, gwas);
#write.csv(tt, "gwas-enrichment.txt", row.names=FALSE, na="");
#print(attr(tt, "overall_enrichment_p"));
#
#tt2 <- cpg_enrichment(cpgs, met_annot);
#write.csv(tt2, "cpg-feature-enrichment.txt", row.names=FALSE, na="");


enrichment <- function(tot_N, tot_member, num_hits, num_intersect) {
	# Construct contingency matrix and divide them into hits vs. non-hits and member vs. non-members
	# First row: Non-members, Second row: Members.
	# First column: Non-hits, Second column: hits.
	cont_mat <- matrix(c(tot_N - tot_member - num_hits + num_intersect, tot_member - num_intersect, num_hits - num_intersect, num_intersect), nrow=2);
	return(fisher.test(cont_mat, alternative="greater")$p.value);
}

depletion <- function(tot_N, tot_member, num_hits, num_intersect) {
	# Construct contingency matrix and divide them into hits vs. non-hits and member vs. non-members
	# First row: Non-members, Second row: Members.
	# First column: Non-hits, Second column: hits.
	cont_mat <- matrix(c(tot_N - tot_member - num_hits + num_intersect, tot_member - num_intersect, num_hits - num_intersect, num_intersect), nrow=2);
	return(fisher.test(cont_mat, alternative="less")$p.value);
}

# Tally CpG features. Note: These are available in Illumina's standard annotation set for EPIC chips
# tt is methylation annotation table (full or a subset)
tally_cpg_features <- function(tt) {
	tbl <- table(tt[, "Relation_to_Island"]);
	dmr <- table(tt[, "DMR"]);
	n <- c(tbl["Island"], tbl["N_Shelf"], tbl["N_Shore"], tbl["OpenSea"], tbl["S_Shelf"], tbl["S_Shore"],
		length(grep("Body", tt[, "UCSC_RefGene_Group"])),
		length(grep("1stExon", tt[, "UCSC_RefGene_Group"])),
		length(grep("ExonBnd", tt[, "UCSC_RefGene_Group"])),
		length(grep("3'UTR", tt[, "UCSC_RefGene_Group"])),
		length(grep("5'UTR", tt[, "UCSC_RefGene_Group"])),
		length(grep("TSS200", tt[, "UCSC_RefGene_Group"])),
		length(grep("TSS1500", tt[, "UCSC_RefGene_Group"])),
		sum(dmr[names(dmr) != ""]),
		dmr["CDMR"],
		dmr["RDMR"],
		dmr["DMR"]);
	## do DHS
	n <- c(n, sum(tt[, "DNase_Hypersensitivity_NAME"] != ""));
	n <- c(n, sum(tt[, "Phantom5_Enhancers"] != ""));
	n <- c(n, sum(tt[, "TFBS_NAME"] != ""));
	n <- c(n, sum(tt[, "OpenChromatin_NAME"] != ""));	
	 
	names(n) <- c("Island", "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore", "Body", "1stExon", "ExonBnd","3UTR", "5UTR", "TSS200", "TSS1500",
		"DMR_All", "CDMR", "RDMR", "DMR_Other", "DHS", "Phantom5_Enhancer","TFBS","OpenChromatin");
	return(n);
}

# Pre-process standard methylation annotation
# Example usage: library(data.table); met_annot <- prep_annotation(data.frame(fread("F:/Annotation/met-annot.txt")));
prep_annotation <- function(met_annot, gene_name_col = "UCSC_RefGene_Name") {
	attr(met_annot, "tally_features") <- tally_cpg_features(met_annot);
	GeneSymbol <- unlist(lapply(met_annot[,gene_name_col], function(x) { paste(sort(unique(unlist(strsplit(x, ";")))), collapse=";"); }));
	met_annot[, "GeneSymbol"] <- GeneSymbol;
	#met_annot <- met_annot[met_annot[, "chr"] != "", ]; # We don't want test probes - none in annotation for EPIC array.
	return (met_annot);
}

# Prepare the NHGRI GWAS catalog.
# Example usage:
# gwascat <- prep_gwas(read.delim("~/projects/nhgri-gwas/gwascatalog-20150326.txt", header=TRUE, sep="\t", as.is=TRUE));
prep_gwas <- function(gwascat, t0=5e-8) {
	# New versions of GWAS catalog seems to sport upper case column names.
	# Make them uniform. Note, though, that the new version of p.Value mLog is PVALUE_MLOG, which I don't use here anyway.
	colnames(gwascat) <- toupper(colnames(gwascat));
	snp_colpos <- which(colnames(gwascat) == "SNPS");
	pval_colpos <- which(colnames(gwascat) == "P.VALUE");
	all_gwas_snps <- strsplit(gwascat[,snp_colpos], ", ");
	num_snps <- unlist(lapply(all_gwas_snps, length));
	tbl <- c();
	for (i in (1:NCOL(gwascat))[-snp_colpos]) tbl <- cbind(tbl, rep(gwascat[,i], num_snps));
	all_gwas_snps <- unlist(all_gwas_snps);
	tbl <- cbind(all_gwas_snps, tbl);
	colnames(tbl) <- c("SNPS", colnames(gwascat)[-snp_colpos]);
	tbl[tbl[,pval_colpos] == "E", pval_colpos] <- "1E-6";
	tbl[tbl[,pval_colpos] %in% c("NS", "Pending"), pval_colpos] <- NA;
	tbl[,pval_colpos] <- as.numeric(tbl[,pval_colpos]);
	tbl <- tbl[!is.na(tbl[, pval_colpos]), ];

	x <- tbl[, "CHR_ID"];
	x[x == "NR"] <- NA;
	x[x == "23"] <- "X";
	tbl[, "CHR_ID"] <- x;
	x <- tbl[, "CHR_POS"];
	x[x == "NR"] <- NA;
	tbl[, "CHR_POS"] <- x;

	# Some p-values are mangled (mixed up signs). Fixed at November 3, 2015
	pval <- as.numeric(tbl[, pval_colpos]);
	pval[pval > 1] <- 10^(-log10(pval[pval > 1]));
	tbl[, pval_colpos] <- pval;

	tbl <- tbl[pval <= t0, ];

	return(tbl);
}

# CpG feature enrichment.
# cpg_list is a vector of significant CpG results.
# met_annot is the pre-treated methylation annotation (using prep_annotation)
cpg_enrichment <- function (cpg_list, met_annot) {
	cpg_tot <- NROW(met_annot);
	tbl <- met_annot[met_annot[, "Name"] %in% cpg_list, ];
	cpg_cur <- NROW(tbl);
	stopifnot(cpg_cur > 0);

	if (is.null(attr(met_annot, "tally_features"))) {
		attr(met_annot, "tally_features") <- tally_cpg_features(met_annot);
	}
	n_all <- attr(met_annot, "tally_features");
	n_hit <- tally_cpg_features(tbl);
	ee <- rep(1, length(n_all));
	dd <- rep(1, length(n_all));
	for (i in 1:length(n_all)) {
		#cat(i,":", names(n_all)[i], ":", n_hit[i], "\n");
		if (is.na(n_hit[i])) n_hit[i] = 0;
		ee[i] <- enrichment(cpg_tot, n_all[i], cpg_cur, n_hit[i]);
		dd[i] <- depletion(cpg_tot, n_all[i], cpg_cur, n_hit[i]);
	}
	cpg_tbl <- data.frame(Enrichment_P=ee, Depletion_P=dd);
	rownames(cpg_tbl) <- names(n_all);
	return(cpg_tbl);
}

# GWAS enrichment analysis
# cpg_list is a vector of significant CpG results.
# met_annot is the pre-treated methylation annotation (using prep_annotation)
# gwascat is the prepared GWAS catalog from the NHGRI (see prep_gwas)
# traits is a vector of traits / diseases you want to perform enrichment on
cpg_gwas_enrichment <- function(cpg_list, met_annot, gwascat, traits=NULL) {
	if (!("GeneSymbol" %in% colnames(met_annot))) { # Pre-treat if necessary
		GeneSymbol <- unlist(lapply(met_annot[,"UCSC_RefGene_Name"], function(x) { paste(sort(unique(unlist(strsplit(x, ";")))), collapse=";"); }));
		met_annot <- cbind(met_annot, GeneSymbol);
		# Get rid of QC probes, if there are any
		met_annot <- met_annot[met_annot[, grep("Chr", colnames(met_annot), TRUE)] != "",];
		rm(GeneSymbol);
	}
	tbl <- met_annot[met_annot[, "Name"] %in% cpg_list, ];
	rep_gene_col <- which(colnames(gwascat) == "REPORTED.GENE.S.");
	map_gene_col <- which(colnames(gwascat) == "MAPPED_GENE");
	disease_col <- which(colnames(gwascat) == "DISEASE.TRAIT");
	genesym_col <- which(colnames(met_annot) == "GeneSymbol");
	
	if (!is.null(traits)) {
		gwascat <- gwascat[gwascat[, disease_col] %in% traits, ];
	}
	stopifnot(NROW(gwascat) > 0);

	met_genes <- unique(unlist(strsplit(met_annot[, genesym_col], ";")));
	num_unique_genes <- length(met_genes);
	# 27322

	# Find GWAS genes that are only in the platform
	gene_gwas1 <- sort(gsub(" ", "", unique(unlist(strsplit(gwascat[, rep_gene_col], ",")))));
	gene_gwas2 <- sort(gsub(" ", "", unique(unlist(strsplit(gwascat[, map_gene_col], " - ")))));
	gene_gwas <- sort(intersect(union(gene_gwas1, gene_gwas2), met_genes));
	N_gwas <- length(gene_gwas);
	b_gwas <- unlist(lapply(tbl[, genesym_col], function(x) sum(unique(unlist(strsplit(x, ";"))) %in% gene_gwas) > 0));

	gene_meta <- unique(unlist(strsplit(tbl[, genesym_col], ";")));
	b_meta1 <- unlist(lapply(gwascat[, rep_gene_col], function(x) sum(unique(unlist(strsplit(x, ", "))) %in% gene_meta) > 0));
	b_meta2 <- unlist(lapply(gwascat[, map_gene_col], function(x) sum(unique(unlist(strsplit(x, " - "))) %in% gene_meta) > 0));
	b_meta <- b_meta1 | b_meta2;

	hit_genes <- unique(unlist(strsplit(tbl[, genesym_col], ";")));
	N_hit <- length(hit_genes);
	N_intersect <- length(unique(tbl[b_gwas, genesym_col]));
	overall_p <- enrichment(num_unique_genes, N_gwas, N_hit, N_intersect);
	diseases <- sort(unique(gwascat[b_meta, disease_col]));
	n_genes_hit <- rep(0, length(diseases));
	n_genes_in_gwas <- rep(0, length(diseases));
	pval <- rep(0, length(diseases));
	genes_in_gwas <- rep("", length(diseases));
	genes_in_list <- rep("", length(diseases));
	
	for (i in 1:length(diseases)) {
		disease <- diseases[i];
		b <- gwascat[, disease_col] == disease;
		gene_gwas1_d <- sort(gsub(" ", "", unique(unlist(strsplit(gwascat[b, rep_gene_col], ",")))));
		gene_gwas2_d <- sort(gsub(" ", "", unique(unlist(strsplit(gwascat[b, map_gene_col], " - ")))));
		gene_gwas_d <- sort(unique(intersect(union(gene_gwas1_d, gene_gwas2_d), met_genes)));
		n1 <- length(gene_gwas_d);

		#b2 <- unlist(lapply(tbl[, genesym_col], function(x) sum(unique(unlist(strsplit(x, ";"))) %in% gene_gwas_d) > 0));
		#n2 <- sum(b2);
		#n2_alt <- length(unique(intersect(gene_gwas_d, gene_meta)));
		#if (n2 > n2_alt) n2 <- n2_alt;
		#if (n2 > n1) n2 <- n1;
	
		hits <- unique(intersect(gene_gwas_d, gene_meta));
		n2 <- length(hits);
		n_genes_in_gwas[i] <- n1;
		n_genes_hit[i] <- n2;
		#cat(i, disease, num_unique_genes, n1, N_hit, n2, "\n");
		#cat(gene_gwas_d, "\n");
		#cat(unique(intersect(gene_gwas_d, gene_meta)), "\n");
		pval[i] <- enrichment(num_unique_genes, n1, N_hit, n2);
		#cat(pval[i], "\n");
		genes_in_gwas[i] <- paste(gene_gwas_d, collapse=";");
		genes_in_list[i] <- paste(hits, collapse=";");
	}
	fdr = p.adjust(pval, "BH");

	gwas_enrichment_tbl <- data.frame(Diseases=diseases, N_Genes_In_GWAS = n_genes_in_gwas, N_Genes_In_Dataset = n_genes_hit, Enrichment_P = pval,
		Enrichment_FDR = fdr, GenesInGWASCat = genes_in_gwas, GenesInList = genes_in_list);
	gwas_enrichment_tbl <- gwas_enrichment_tbl[order(gwas_enrichment_tbl[, "Enrichment_P"]), ];
	attr(gwas_enrichment_tbl, "overall_enrichment_p") <- overall_p;
	attr(gwas_enrichment_tbl, "n_genes_in_gwas") <- N_gwas;
	attr(gwas_enrichment_tbl, "n_genes_hit") <- N_hit;
	attr(gwas_enrichment_tbl, "n_traits") <- length(diseases);
	
	return(gwas_enrichment_tbl);
}
