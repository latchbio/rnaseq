#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("One argument must be supplied (GTF file).n", call. = FALSE)
}

if (!require("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
}
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")

library(tximport)
library(GenomicFeatures)

quant_path <- args[2]

# NEEDS TO BE A ARGUMENT
gtf_path <- args[3]

txdb <- makeTxDbFromGFF(gtf_path)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi <- tximport(quant_path, type = "salmon", tx2gene = tx2gene)

counts <- txi$counts

df <- data.frame(Name = row.names(counts), NumReads = counts)

write.table(df, file = args[4], quote = FALSE, row.names = FALSE, sep = "\t")
