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
labels <- gsub(".normalised.*$", "",
               names(data)[grepl(".normalised.*$", names(data))])
colours <- as.numeric(samples$condition)
for ( cluster in unique(sort(data$cluster)) ) {
    data.subset <- data[data$cluster == cluster,]
    write.table(data.subset, file=paste0(outputBase, '-cluster-', cluster,
                                         '.tsv'),
                quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    # Plot counts
    pdf(paste0(outputBase, '-cluster-', cluster,
               '-counts.pdf'))
    for (i in 1:nrow(data.subset)) {
        counts <- data.subset[i, grepl(".normalised.*$", names(data.subset)) ]
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
        axis(1, at=1:length(labels), lab=labels, las=2, cex.axis=0.5)
        axis(2)
        title(main=sprintf("%s:%d-%d\n%s / %s\n%.2f",
                           data.subset[i,"Chr"], data.subset[i,"Start"],
                           data.subset[i,"End"], data.subset[i,"Gene.ID"],
                           data.subset[i,"Name"], data.subset[i,"adjp"]))
        title(xlab="")
        title(ylab="Normalised Counts")
        legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
               pt.bg=1:length(levels(samples$condition)))
    }
    graphics.off()
}
