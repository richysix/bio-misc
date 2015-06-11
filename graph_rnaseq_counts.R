#!/usr/bin/env Rscript

# Graph count data from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",            Args[6])
samplesFile <- ifelse(is.na(Args[7]), "deseq2/samples.txt", Args[7])
pdfFile     <- ifelse(is.na(Args[8]), "counts.pdf",         Args[8])

# Read data
data <- read.delim(dataFile, header=TRUE)

# Support different column names
names(data)[names(data) == 'chr']     <- 'Chr'
names(data)[names(data) == 'start']   <- 'Start'
names(data)[names(data) == 'end']     <- 'End'
names(data)[names(data) == 'ID']      <- 'Gene.ID'
names(data)[names(data) == 'adjpval'] <- 'adjp'

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Graph parameters
labels <- gsub(".normalised.*$", "",
               names(data)[grepl(".normalised.*$", names(data))])
colours <- as.numeric(samples$condition)

pdf(pdfFile)

# Plot each region separately
for (i in 1:nrow(data)) {
    counts <- data[i, grepl(".normalised.*$", names(data)) ]
    par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
    axis(1, at=1:length(labels), lab=labels, las=2, cex.axis=0.5)
    axis(2)
    title(main=sprintf("%s:%d-%d\n%s\n%.2f",
                       data[i,"Chr"], data[i,"Start"],
                       data[i,"End"], data[i,"Gene.ID"],
                       data[i,"adjp"]))
    title(xlab="")
    title(ylab="Normalised Counts")
    legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
           pt.bg=1:length(levels(samples$condition)))
}

graphics.off()
