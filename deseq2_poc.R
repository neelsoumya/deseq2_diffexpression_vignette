#################################################################
# Differential expression testing using DESeq2 
#
# Usage: nohup R --no-save < deseq2_poc.R
# OR
# R # then
# source("deseq2_poc.R")
#
# Installation:
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("airway")
#
# Adapted from
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#################################################################


##########################################
# Load libraries
##########################################
library(DESeq2)
library("airway")


##########################################
# Load example data
##########################################
data("airway")
se <- airway

##########################################
# The constructor function below shows 
# the generation of a DESeqDataSet from 
# a RangedSummarizedExperiment se.
##########################################
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE


##########################################
# Pre-filtering
# perform a minimal pre-filtering to keep
# only rows that have at least 10 reads 
# total. Note that more strict filtering
# to increase power is automatically
# applied via independent filtering on
# the mean of normalized counts within
# the results function.
##########################################
i_num_reads_per_gene = 10
keep <- rowSums(counts(ddsSE)) >= i_num_reads_per_gene
ddsSE_filtered <- ddsSE[keep,]


##########################################
# Factor levels
##########################################
ddsSE_filtered$conditions = factor(ddsSE_filtered$conditions, 
                                   levels = c("untreated","treated"))


##########################################
# Perform differential abundance test
##########################################
ddsSE_filtered <- DESeq(ddsSE_filtered)
res <- results(ddsSE_filtered)

##########################################
# Output results
##########################################
res
resultsNames(ddsSE_filtered)

lfc_results <- lfcShrink(ddsSE_filtered, coef = 2)
lfc_results

# order results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

# summarize some basic tallies using the summary function.
summary(res)

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(ddsSE_filtered, alpha=0.05)
summary(res05)

# How many adjusted p-values were less than 0.05?
sum(res05$padj < 0.05, na.rm=TRUE)


##########################################
# Visualize results
# Points will be colored red if the adjusted p value is less than 0.1.
# Points which fall out of the window are plotted as open triangles pointing either up or down.
##########################################
plotMA(res)
plotMA(res$log2FoldChange)

# plot counts
plotCounts(ddsSE_filtered, gene=which.min(res$padj), intgroup="condition")


##########################################
# Save results
##########################################
write.csv(as.data.frame(resOrdered), 
          file="ordered_results_by_condition.csv")



##########################################
# Miscellaneous
##########################################

# Multi-factorial design
colData(ddsSE_filtered)

ddsSE_filtered_MF <- ddsSE_filtered
levels(ddsSE_filtered_MF$type)

levels(ddsSE_filtered_MF$type) <- sub("-.*", "", levels(ddsSE_filtered_MF$type))
levels(ddsSE_filtered_MF$type)

# rerun DESEq2
design(ddsSE_filtered_MF) <- formula(~ type + condition)
ddsSE_filtered_MF <- DESeq(ddsSE_filtered_MF)

