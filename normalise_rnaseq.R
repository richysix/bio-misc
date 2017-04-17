#!/usr/bin/env Rscript

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))

Args        <- commandArgs()
countFile   <- ifelse(is.na(Args[6]),  "deseq2/counts.txt",  Args[6])
samplesFile <- ifelse(is.na(Args[7]),  "deseq2/samples.txt", Args[7])
outputDir   <- ifelse(is.na(Args[8]),  "deseq2",             Args[8])

# Get data and samples
countData <- read.table(   countFile, header=TRUE, row.names=1 )
samples   <- read.table( samplesFile, header=TRUE, row.names=1 )

# Normalise
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ 1)
dds <- estimateSizeFactors(dds)

# Write out normalisation size factors
write.table(sizeFactors(dds),
            file=paste(outputDir, "/size-factors.txt", sep=""),
            col.names=FALSE, quote=FALSE, sep="\t")

# Write out normalised counts
write.table(counts(dds, normalized=TRUE),
            file=paste(outputDir, "/normalised-counts.txt", sep=""),
            col.names=FALSE, quote=FALSE, sep="\t")
