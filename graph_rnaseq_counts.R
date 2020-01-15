#!/usr/bin/env Rscript

# Graph count data from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

library(ggplot2)
library(reshape2)
suppressPackageStartupMessages(library(dplyr))
library(svglite)

Args           <- commandArgs()
dataFile       <- ifelse(is.na(Args[6]),  "all.tsv",            Args[6])
samplesFile    <- ifelse(is.na(Args[7]),  "deseq2/samples.txt", Args[7])
outputFile     <- ifelse(is.na(Args[8]),  "counts.pdf",         Args[8])
plotStyle      <- ifelse(is.na(Args[9]),  "default",            Args[9])
shapeVariable  <- ifelse(is.na(Args[10]), "none",               Args[10])
colourVariable <- ifelse(is.na(Args[11]), "stage",              Args[11])

# Read data
data <- read.delim(dataFile, header=TRUE, check.names=FALSE)

# Support different column names
names(data)[names(data) == 'chr']     <- 'Chr'
names(data)[names(data) == 'start']   <- 'Start'
names(data)[names(data) == 'end']     <- 'End'
names(data)[names(data) == 'ID']      <- 'Gene ID'
names(data)[names(data) == 'adjpval'] <- 'adjp'

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Check shape and colour variables exist in the samples file
if (shapeVariable != 'none') {
    if (!any(grepl(shapeVariable, colnames(samples)))) {
        stop(paste0('The shape variable, ', shapeVariable,
                    ' does not exist as a column in the samples file'))
    }
    # Check colourVariable as well
    if (!any(grepl(colourVariable, colnames(samples)))) {
        stop(paste0('The colour variable, ', colourVariable,
                    ' does not exist as a column in the samples file'))
    }
}

# Get counts
countData <- data[,grepl(" normalised.*$", names(data))]
names(countData) <- gsub(" normalised.*$", "", names(countData))

# Subset and reorder count data
countData <- countData[, row.names(samples)]

# Graph parameters
shapePalette <- 21:25
if (shapeVariable == 'none') { # Don't use shape
    colours <- as.numeric(samples$condition)
    shapes <- rep(shapePalette[1], length(samples$condition))
} else {
    # Check there aren't too many levels of condition for shape
    if (nlevels(samples[[shapeVariable]]) > length(shapePalette)) {
        stop(paste0('The shape variable, ', shapeVariable,
                    ' has more levels than available shapes (',
                    length(samples$condition), ')'))
    } else {
        shapes <- shapePalette[ as.numeric(samples[[shapeVariable]]) ]
        colours <- as.numeric(samples[[colourVariable]])
    }
}

if (grepl("pdf$", outputFile)) {
    pdf(outputFile)
}

if (grepl("violin", plotStyle)) {
    countData$id <- row.names(countData)
    counts <- melt(countData, id.vars="id", variable.name="condition",
                   value.name="count")
    countData$name <- sprintf("%s:%d-%d\n%s / %s\n%.2f",
                              data[,"Chr"],
                              data[,"Start"],
                              data[,"End"],
                              data[,"Gene ID"],
                              data[,"Name"],
                              data[,"adjp"])
    levels(counts$condition) <- samples$condition
    for (i in 1:nrow(data)) {
        p <- ggplot(counts[counts$id == i,],
               aes(x=condition, y=count, color=condition)) +
            geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) +
            labs(x="", y="Normalised Counts", title=countData$name[i]) +
            theme_bw() +
            theme(legend.position='none',
                  plot.title = element_text(hjust = 0.5, face="bold"))
        if (plotStyle == "violinpoints") {
            p <- p + geom_jitter()
        }
        if (grepl("pdf$", outputFile)) {
            print(p)
        } else {
            ggsave(filename=paste0(outputFile, i, '.svg'), plot=p, width=10,
                   height=8)
        }
    }
} else {
    for (i in 1:nrow(data)) {
        if (!grepl("pdf$", outputFile)) {
            if (grepl("eps$", outputFile)) {
                epsFile <- gsub('eps$', paste0(i, '.eps'), outputFile)
                postscript(file = epsFile, width = 7, height = 6,
                           paper = 'special', horizontal = FALSE)
            } else {
                # Assume SVG
                svglite(paste0(outputFile, i, '.svg'))
            }
        }
        if (shapeVariable == 'none') {
            palette(rainbow(length(levels(samples$condition))))
        } else {
            palette(rainbow(length(levels(samples[[colourVariable]]))))
        }
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(countData[i,]), axes=FALSE, ann=FALSE, pch=shapes,
             bg=colours)
        axis(1, at=1:ncol(countData), lab=names(countData), las=2, cex.axis=0.5)
        axis(2)
        title(main=sprintf("%s:%d-%d\n%s / %s\n%.2f",
                           data[i,"Chr"], data[i,"Start"],
                           data[i,"End"], data[i,"Gene ID"],
                           data[i,"Name"], data[i,"adjp"]))
        title(xlab="")
        title(ylab="Normalised Counts")

        if (shapeVariable == 'none') {
            legend("topright", inset=c(0, -0.1), levels(samples$condition),
                   pch=21, pt.bg=1:length(levels(samples$condition)))
        } else {
            legend("topright", inset=c(0, -0.1),
                   levels(samples[[shapeVariable]]),
                   pch=shapePalette[ seq_len(nlevels(samples[[shapeVariable]])) ])
            legend("right", inset=c(0, -0.1), levels(samples[[colourVariable]]),
                   pch=21, pt.bg=1:length(levels(samples[[colourVariable]])))
        }
        if (!grepl("pdf$", outputFile)) {
            graphics.off()
        }
    }
}

if (grepl("pdf$", outputFile)) {
    graphics.off()
}
