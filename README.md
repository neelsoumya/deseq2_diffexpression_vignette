# deseq2_diffexpression_vignette

Differential expression analysis of single-cell sequencing data using DESeq2.

 	Adapted from
	http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


* Installation:

	In R
	
	NOTE: try http:// if https:// URLs are not supported
	
	
	source("https://bioconductor.org/biocLite.R")
	

	biocLite("airway")
	
	if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
     BiocManager::install("DESeq2")


* Usage:


	nohup R --no-save < deseq2_poc.R
	

	OR

	in R 

	source("deseq2_poc.R")



