#!/usr/bin/env Rscript

# Hard clustering and cluster size estimation from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

suppressPackageStartupMessages(library(hopach))
library(naturalsort)

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",            Args[6])
samplesFile <- ifelse(is.na(Args[7]), "deseq2/samples.txt", Args[7])
outputBase  <- ifelse(is.na(Args[8]), "all",                Args[8])

# Read data
data <- read.delim(dataFile, header=TRUE, check.names=FALSE)

# Support different column names
names(data)[names(data) == 'chr']     <- 'Chr'
names(data)[names(data) == 'start']   <- 'Start'
names(data)[names(data) == 'end']     <- 'End'
names(data)[names(data) == 'ID']      <- 'Gene.ID'
names(data)[names(data) == 'adjpval'] <- 'adjp'

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Ensure chromosome is a factor, even if all numeric
data[,1] <- factor(data[,1])

# Get counts
countData <- data[,grepl(".normalised.counts?$", names(data))]
names(countData) <- gsub(".normalised.counts?$", "", names(countData))

# Get median counts
medianData <- matrix(nrow=nrow(countData),
                     ncol=length(levels(samples$condition)),
                     dimnames=list(rownames(countData),
                                   naturalsort(levels(samples$condition))))
for (condition in naturalsort(levels(samples$condition))) {
    medianData[,condition] <-
        apply(countData[,samples$condition == condition, drop=FALSE], 1,
              median)
}

# Cluster
distMatrix <- distancematrix(medianData, 'cosangle')
hobj <- hopach(medianData, dmat=distMatrix)

# Output clusters and number of clusters
write.table(hobj$clust$k, file=paste0(outputBase, '-num-clusters.tsv'),
            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
options(scipen=100)
clusters <- as.character(hobj$clust$labels)
options(scipen=0)
data$cluster <- clusters
write.table(data, file=paste0(outputBase, '-clusters.tsv'),
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
