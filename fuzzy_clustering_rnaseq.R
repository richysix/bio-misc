#!/usr/bin/env Rscript

# Fuzzy clustering from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(Mfuzz))

Args           <- commandArgs()
dataFile       <- ifelse(is.na(Args[6]),  "all.tsv",            Args[6])
samplesFile    <- ifelse(is.na(Args[7]),  "deseq2/samples.txt", Args[7])
outputBase     <- ifelse(is.na(Args[8]),  "all",                Args[8])
numClusters    <- ifelse(is.na(Args[9]),  100,       as.integer(Args[9]))
alphaThreshold <- ifelse(is.na(Args[10]), 0.7,       as.numeric(Args[10]))

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
                                   levels(samples$condition)))
for (condition in levels(samples$condition)) {
    medianData[,condition] <- apply(countData[,samples$condition == condition],
                                    1, median)
}

# Standardise data
eset <- ExpressionSet(assayData=as.matrix(medianData))
eset <- filter.std(eset, min.std=0, visu=FALSE)
eset.s <- standardise(eset)

# Get fuzzy clusters
pdf(paste0(outputBase, '-', numClusters, '-', alphaThreshold, '-mfuzz.pdf'))
cl <- mfuzz(eset.s, c=numClusters, m=mestimate(eset.s))
mfuzz.plot(eset.s, cl=cl, mfrow=c(4,4), new.window=FALSE,
           time.labels=colnames(medianData))
graphics.off()

# Output fuzzy clusters
acore.list <- acore(eset.s, cl=cl, min.acore=alphaThreshold)
for (i in 1:numClusters) {
    data.subset <- data[c(as.character(acore.list[[i]][,1])),]
    write.table(data.subset, file=paste0(outputBase, '-', numClusters, '-',
                                         alphaThreshold, '-', i, '.tsv'),
                quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
}
