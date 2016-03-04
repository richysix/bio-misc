#!/usr/bin/env Rscript

# Volcano plot from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",     Args[6])
outputFile  <- ifelse(is.na(Args[7]), "volcano.png", Args[7])

# Read data
data <- read.delim(dataFile, header=TRUE, check.names=FALSE)

# Support different column names
names(data)[names(data) == 'adjpval'] <- 'adjp'

# Plot (red if adjp < 0.05; orange if log2fc > 1; green if both)
png(outputFile)
with(data, plot(log2fc, -log10(adjp), pch=20, main="Volcano plot",
                xlim=c(-max(abs(log2fc)), max(abs(log2fc)))))
with(subset(data, adjp < 0.05 ), points(log2fc, -log10(adjp), pch=20,
                                        col="red"))
with(subset(data, abs(log2fc) > 1), points(log2fc, -log10(adjp), pch=20,
                                           col="orange"))
with(subset(data, adjp < 0.05 & abs(log2fc) > 1), points(log2fc, -log10(adjp),
                                                         pch=20, col="green"))
graphics.off()
