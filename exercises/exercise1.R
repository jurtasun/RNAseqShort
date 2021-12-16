# RNAseq course - Exercise 1
# Jesus Urtasun - December 2021
print("RNAseq course - Exercise 1")
print("Jesus Urtasun - December 2021")

# In this exercise, we will read in a count data from erythroblast differentiation experiment in mice.
# Data downloaded from GEO database (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49843) 
# Dara aligned to mm9 genome assembly using Rsubread R package.

# 1
# Read the sample information and count data
# Identify how many factors this experiments involves
setwd("/Volumes/bioinfomatics$/jurtasun/Courses/RNAseqShort/exercises")
suppressPackageStartupMessages(library(DESeq2))

# Load sample description generated for this exercise
targetsE <- read.table("Exercise_ShortRNAseq_sample.info", sep = "\t", header = TRUE)
# targetsE

# Load the count information
all_countsE <- read.csv(file = "Exercise_ShortRNAseq_counts.csv", header = T, row.names = 1)
# head(all_countsE)

# In this experiment we have two factors: condition and batch
# For the condition factor there are 3 levels: FFa, KOa and KOb
# For the batch factor there are 3 level: a, b and c

# 2
# Construct DESEq dataset object using sampple information and counts data
# Provide the entrez_id as identifier for the exercises
cDataE <- data.frame(name = targetsE$sample, condition = targetsE$condition, batch = targetsE$batch)
ddsE <- DESeqDataSetFromMatrix(countData = all_countsE, colData = cDataE, design = ~condition)

# 3
# Find the number of genes that are changed in knockdown samples versus control, irrespective of fold change
# Find the number of genes that are changed in the above situation with fold change threshold, fold change ratio > 2
ddsE <- DESeq(ddsE)

# Show how many genes were differentially expressed in KOa vs FFa, with FDR < 0.05
res1 <- results(ddsE, contrast = c("condition", "KOa", "FFa"))
# > summary(res1, alpha = 0.05)

# Shows how many genes were differentially expressed in KOa vs FFa, whith FDR < 0.05 and fold change ratio > 2
DE_res1 <- res1[complete.cases(res1$padj), ]
DE_res1 <- DE_res1[DE_res1$padj <- 0.05 & abs(DE_res1$log2FoldChange) > 1, ]

# Show how many genes were differentially expressed in KOa vs FFa, with FDR < 0.05
res2 <- results(ddsE, contrast = c("condition", "KOb", "FFa"))
# > summary(res2, alpha = 0.05)

# Show how many genes were differentially expressed in KOa vs FFa, with FDR < 0.05 and fold change ratio > 2
DE_res2 <- res2[complete.cases(res2$padj), ]
DE_res2 <- DE_res2[DE_res2$padj <- 0.05 & abs(DE_res2$log2FoldChange) > 1, ]

# 4
# Using Biomart - Add extra columns to the results dataframe such as gene name, gene biotype, gene description
# library(biomaRt)
# 
# mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2012.archive.ensembl.org")
# 
# bm <- getBM(attributes = c("ensembl_gene_id", "gene_biotype" ,"description" ,"mgi_symbol"),
#             filters = "entrezgene", values = rownames(res1), mart=mart)
# 
# resAnnotated <- merge(as.data.frame(res1),bm,by.x=0,by.y=1)