#!/usr/bin/env Rscript

# PCA from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

#suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

Args               <- commandArgs()
dataFile           <- ifelse(is.na(Args[6]),  "all.tsv",            Args[6])
samplesFile        <- ifelse(is.na(Args[7]),  "deseq2/samples.txt", Args[7])
outputFile         <- ifelse(is.na(Args[8]),  "pca.pdf",            Args[8])
outputBase         <- ifelse(is.na(Args[9]),  "all",                Args[9])
transformMethod    <- ifelse(is.na(Args[10]), "rlog",               Args[10])
regionCount        <- ifelse(is.na(Args[11]), 500,         as.integer(Args[11]))
varPCThreshold     <- ifelse(is.na(Args[12]), 1,           as.numeric(Args[12]))
varRegionThreshold <- ifelse(is.na(Args[13]), 0.1,         as.numeric(Args[13]))
includeLoadings    <- ifelse(is.na(Args[14]), FALSE,       as.logical(Args[14]))

remove_common_prefix <- function(names) {
    common_prefix_length <- 0
    got_longest <- FALSE
    while ( !got_longest ) {
        prefixes <- substr(names, 1, common_prefix_length + 1)
        if ( length(unique(prefixes)) == 1 ) {
            common_prefix_length <- common_prefix_length + 1
        } else {
            got_longest <- TRUE
        }
    }
    common_prefix <- substr(names[1], 1, common_prefix_length)
    names <- gsub(common_prefix, "", names)
    return(names)
}

# Read data
data <- read.delim(dataFile, header=TRUE, check.names=FALSE)

# Support different column names
names(data)[names(data) == 'chr']     <- 'Chr'
names(data)[names(data) == 'start']   <- 'Start'
names(data)[names(data) == 'end']     <- 'End'
names(data)[names(data) == 'ID']      <- 'Gene.ID'
names(data)[names(data) == 'adjpval'] <- 'adjp'

# Read samples
samples <- read.delim( samplesFile, header=TRUE, row.names=1 )

# Ensure chromosome is a factor, even if all numeric
data[,'Chr'] <- factor(data[,'Chr'])

# Get counts
countData <- data[,grepl(".count$", names(data)) &
                  !grepl(".normalised.count$", names(data))]
names(countData) <- gsub(".count$", "", names(countData))

# Subset and reorder count data
countData <- countData[, row.names(samples)]

# Transform using DESeq2
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ 1)
dds <- estimateSizeFactors(dds)
if (transformMethod == "rlog") {
    dds <- rlogTransformation(dds, blind=TRUE)
} else {
    dds <- varianceStabilizingTransformation(dds, blind=TRUE)
}

# PCA
rv <- rowVars(assay(dds))
select <- order(rv, decreasing=TRUE)[seq_len(min(regionCount, length(rv)))]
pca <- prcomp(t(assay(dds)[select,]))
propVarPC <- pca$sdev^2 / sum( pca$sdev^2 )
aload <- abs(pca$rotation)
propVarRegion <- sweep(aload, 2, colSums(aload), "/")
if (includeLoadings) {
    loadings <- sweep(pca$rotation, 2, colSums(aload), "/")
}

# Output PC coordinates
lastSigPC <- sum(propVarPC * 100 >= varPCThreshold)
x <- pca$x[,1:lastSigPC]
if (grepl("csv$", dataFile)) {
    write.csv( x, file=paste0(outputBase, "-PCs.csv"))
} else {
    write.table( x, file=paste0(outputBase, "-PCs.tsv"), quote=FALSE, sep="\t" )
}

# Output regions contributing most to each PC
for (i in seq.int(lastSigPC)) {
    data[select, "% variance explained"] <- propVarRegion[,i] * 100
    if (includeLoadings) {
        data[select, "loadings"] <- loadings[,i]
    }
    topData <- subset(data[order(data$`% variance explained`,
        decreasing=TRUE),], `% variance explained` >= varRegionThreshold)
    write.table( topData, file=paste0(outputBase, "-PC", i, ".tsv"),
        row.names=FALSE, quote=FALSE, sep="\t" )
}

# Write PCs to PDF
if (grepl("pdf$", outputFile)) {
    pdf(outputFile)
}

# Show variance explained for each PC
d <- data.frame(PC=seq.int(length(propVarPC)), var=propVarPC)
p <- ggplot(data=d, aes(x=PC, y=var)) + geom_line() + geom_point() +
    xlab("PC") + ylab("Variance explained") + ylim(c(0, 1))
if (grepl("pdf$", outputFile)) {
    print(p)
} else {
    ggsave(filename=paste0(outputFile, 'expl.svg'), plot=p, width=10,
           height=8)
}
p <- ggplot(data=d, aes(x=PC, y=cumsum(var))) + geom_line() + geom_point() +
    xlab("PC") + ylab("Cumulative variance explained") + ylim(c(0, 1))
if (grepl("pdf$", outputFile)) {
    print(p)
} else {
    ggsave(filename=paste0(outputFile, 'cum.expl.svg'), plot=p, width=10,
           height=8)
}

# Plot PCs in pairs
intgroup.df <- as.data.frame(colData(dds)[, c("condition"), drop=FALSE])
group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))
short_sample_names <- remove_common_prefix(colnames(dds))
for (i in seq.int(lastSigPC - 1)) {
    first <- i
    second <- i + 1
    d <- data.frame(first=pca$x[,first], second=pca$x[,second], group=group,
        intgroup.df, name=colnames(dds))
    p <- ggplot(data=d, aes_string(x="first", y="second", color="group")) +
        geom_point(size=2) +
        geom_text(aes(label=short_sample_names), hjust=0, vjust=0, size=4,
                  show.legend=FALSE) +
        xlab(paste0("PC", first, ": ", round(propVarPC[first] * 100, 1),
            "% variance")) +
        ylab(paste0("PC", second, ": ", round(propVarPC[second] * 100, 1),
            "% variance"))
    if (grepl("pdf$", outputFile)) {
        print(p)
    } else {
        ggsave(filename=paste0(outputFile, first, '-', second, '.svg'), plot=p,
               width=10, height=8)
    }
}

# Plot proportion of variance explained by each region across chromosome
varMax <- max(propVarRegion[,1:lastSigPC])
propVarRegion <- as.data.frame(propVarRegion[,1:lastSigPC])
propVarRegion$region <- rowMeans(data[select, c('Start', 'End')])
propVarRegion$chr <- data[select, 'Chr']
chrs <- sort(unique(as.numeric(levels(data[,'Chr'])[grepl("^[0-9]+$",
    levels(data[,'Chr']))])))
for (chr in chrs) {
    var_long <- melt(propVarRegion[propVarRegion$chr == chr,1:ncol(propVarRegion)-1],
        id="region", variable.name="PC")
    p <- ggplot(data=var_long, aes(x=region, y=value, colour=PC)) +
        geom_line() + geom_point() +
        xlab(paste0("Chromosome ", chr)) +
        ylab("Variance explained") +
        ylim(c(0, varMax))
    if (grepl("pdf$", outputFile)) {
        print(p)
    } else {
        ggsave(filename=paste0(outputFile, 'chr', chr, '.svg'), plot=p,
               width=10, height=8)
    }
}

if (grepl("pdf$", outputFile)) {
    graphics.off()
}
