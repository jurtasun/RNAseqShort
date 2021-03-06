# RNAseq course - Chapter 1
# Jesus Urtasun - December 2021
print("RNAseq course - Chapter 1")
print("Jesus Urtasun - December 2021")

# Import packaged
library(DESeq2)

# Read the targets file
targets <- read.table("targets.txt", sep = "\t", header = TRUE)

# Load count data
all_counts <- read.csv(file = "AllCounts.csv", row.names = 1)

# Prepare DESeq dataset object - collect sample information
c_data <- data.frame(name = targets$Sample, Group = targets$Group, Batch = targets$Batch)
rownames(c_data) <- c_data[, 1]

# Construct DESeq dataset object - class DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = all_counts, colData = c_data, design =~ Group)

# Standard differential expression analysis steps are wrapped into a single function, DESeq
# DESeq performs normalization, fitting to the model and statistical testing
dds <- DESeq(dds)

# Functions run by DESeq
# estimateSizeFactors()
# estimateDispersions()
# nbinomWaldTest()

# 1. Estimate size factors
# Anders S and Huber W. (2010). Differential expression analysis for sequence count data. Genome Biol.:11(10):R106.
# The sizeFactors vector assigns to each column the count matrix a value, the size factor, such that
# count values in the columns can be brought to a common scale by dividing by the corresponding factor
# > sizeFactors(dds)

# 2. Extimate dispersions
# The function estimateDispersions() obtains gene-wide dispersion estimates
# Then a curve is fit to the estimates to capture the overall trend of dispersion-mean dependence
# > head(dispersions(dds))
# > plotDispEst(dds)

# 3. Binomial Wald test
# Wald test for computing p-values, testing whether each model coefficient differs significantly from zero
# > nbinomWaldTest(dds)
