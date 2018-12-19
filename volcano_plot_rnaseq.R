#!/usr/bin/env Rscript

# Volcano plot from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a
library('svglite')

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",     Args[6])
outputFile  <- ifelse(is.na(Args[7]), "volcano.png", Args[7])
maxLog2Fc   <- suppressWarnings(ifelse(is.na(Args[8]), NA,
                                as.numeric(Args[8])))
maxLogPVal  <- suppressWarnings(ifelse(is.na(Args[9]), NA,
                                as.numeric(Args[9])))
pValLine    <- suppressWarnings(ifelse(is.na(Args[10]), NA,
                                as.numeric(Args[10])))
markedGene  <- ifelse(is.na(Args[11]), NA, Args[11])

# Read data
data <- read.delim(dataFile, header=TRUE, check.names=FALSE)

# Support different column names
names(data)[names(data) == 'adjpval'] <- 'adjp'
names(data)[names(data) == 'Gene ID'] <- 'GeneID'

# Plot (red if adjp < 0.05; orange if log2fc > 1; green if both)
# decide on output type. Default is png
if ( grepl("svg$", outputFile) ) {
    svglite(outputFile)
} else {
    png(outputFile)
}
if ( is.na(maxLog2Fc) ) {
    maxLog2Fc <- c(-max(abs(data$log2fc), na.rm=TRUE),
                   max(abs(data$log2fc), na.rm=TRUE))
} else {
    maxLog2Fc <- c(-1 * maxLog2Fc, maxLog2Fc)
}
if ( is.na(maxLogPVal) ) {
    maxLogPVal <- NULL
} else {
    maxLogPVal <- c(0, maxLogPVal)
}
with(data, plot(log2fc, -log10(adjp), pch=20, main="Volcano plot",
                xlim=maxLog2Fc, ylim=maxLogPVal))
with(subset(data, adjp < 0.05 ), points(log2fc, -log10(adjp), pch=20,
                                        col="red"))
with(subset(data, abs(log2fc) > 1), points(log2fc, -log10(adjp), pch=20,
                                           col="orange"))
with(subset(data, adjp < 0.05 & abs(log2fc) > 1), points(log2fc, -log10(adjp),
                                                         pch=20, col="green"))
if ( !is.na(markedGene) ) {
    with(subset(data, GeneID == markedGene), points(log2fc, -log10(adjp),
                                                    pch=21, col="blue", cex=2))
}
if ( !is.na(pValLine) ) {
    abline(h=-log10(pValLine), col="gray60", lty=2)
}
graphics.off()
