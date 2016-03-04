#!/usr/bin/env Rscript

# Hard clustering and cluster size estimation from RNA-Seq output produced by:
# https://gist.github.com/iansealy/2dca28d07c0764e014df
# or https://gist.github.com/iansealy/b9cbc56bd1affe10d37a

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(hopach))
library(naturalsort)

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",            Args[6])
samplesFile <- ifelse(is.na(Args[7]), "deseq2/samples.txt", Args[7])
outputBase  <- ifelse(is.na(Args[8]), "all",                Args[8])

# Temporarily include modified version of mfuzz.plot2
mfuzz.plot2.tmp <- function(eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels,
    time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes",
    x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", col.lab = "black",
    col.main = "black", col.sub = "black", col = "black", centre = FALSE, centre.col = "black",
    centre.lwd = 2, Xwidth = 5, Xheight = 5, single = FALSE, ...) {
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])

    if (missing(colo)) {
        colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", "#FF6000", "#FF7800",
            "#FF8F00", "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00",
            "#C7FF00", "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70",
            "#00FF87", "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF",
            "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF",
            "#0028FF", "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF",
            "#8000FF", "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048",
            "#FF0030", "#FF0018")
    } else {
        if (colo == "fancy") {
            fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), rep(0, length(c(255:150))))
            fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
            fancy.red <- c(c(0:255), rep(255, length(c(255:0))), c(255:150))
            colo <- rgb(b = fancy.blue/255, g = fancy.green/255, r = fancy.red/255)
        }
    }
    colorseq <- seq(0, 1, length = length(colo))

    for (j in 1:dim(cl[[1]])[[1]]) {
        if (single)
            j <- single
        tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
            if (x11)
                X11(width = Xwidth, height = Xheight)
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            } else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            } else {
                ylim <- ylim.set
            }
            if (!is.na(sum(mfrow))) {
                par(mfrow = mfrow, bg = bg, col.axis = col.axis, col.lab = col.lab,
                  col.main = col.main, col.sub = col.sub, col = col)
            } else {
                par(bg = bg, col.axis = col.axis, col.lab = col.lab, col.main = col.main,
                  col.sub = col.sub, col = col)
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points)))
                xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
        } else {
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            } else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            } else {
                ylim <- ylim.set
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points)))
                xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
        }

        if (length(tmpmem) > 0) {
            for (jj in 1:(length(colorseq) - 1)) {
                tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    if (missing(time.points)) {
                      lines(tmp[tmpind[k], ], col = colo[jj])
                    } else lines(time.points, tmp[tmpind[k], ], col = colo[jj])
                  }
                }
            }
        }
        if (centre) {
            lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
        }
        if (single)
            return()
    }
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
options(scipen=100)
text_clusters <- as.character(hobj$clust$labels)
options(scipen=0)
data$cluster <- text_clusters
clusters <- unique(sort(data$cluster))

# Make data like Mfuzz so as to use Mfuzz's plotting
eset <- ExpressionSet(assayData=as.matrix(medianData))
eset <- filter.std(eset, min.std=0, visu=FALSE)
eset.s <- standardise(eset)
cl <- list()
cl[[1]] <- matrix(, nrow=hobj$clust$k, ncol=ncol(medianData))
cl[[3]] <- as.integer(factor(data$cluster))
cl[[4]] <- matrix(, nrow=nrow(medianData), ncol=0)
for ( cluster in clusters ) {
    cl[[4]] <- cbind(cl[[4]],
                     matrix(as.integer(data$cluster == cluster), ncol=1))
}

# Plot clusters
pdf(paste0(outputBase, '-hopach.pdf'))
mfuzz.plot2.tmp(eset.s, cl=cl, mfrow=c(4,4), x11=FALSE, centre=FALSE,
           time.labels=colnames(medianData), min.mem=1)
graphics.off()

# Output clusters and number of clusters
write.table(hobj$clust$k, file=paste0(outputBase, '-num-clusters.tsv'),
            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(data, file=paste0(outputBase, '-clusters.tsv'),
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
labels <- gsub(".normalised.*$", "",
               names(data)[grepl(".normalised.*$", names(data))])
colours <- as.numeric(samples$condition)
for ( i in seq.int(length(clusters)) ) {
    data.subset <- data[data$cluster == clusters[i],]
    write.table(data.subset, file=paste0(outputBase, '-cluster-', i, '.tsv'),
                quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    # Plot counts
    pdf(paste0(outputBase, '-cluster-', i, '-counts.pdf'))
    for (j in 1:nrow(data.subset)) {
        counts <- data.subset[j, grepl(".normalised.*$", names(data.subset)) ]
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
        axis(1, at=1:length(labels), lab=labels, las=2, cex.axis=0.5)
        axis(2)
        title(main=sprintf("%s:%d-%d\n%s / %s\n%.2f",
                           data.subset[j,"Chr"], data.subset[j,"Start"],
                           data.subset[j,"End"], data.subset[j,"Gene.ID"],
                           data.subset[j,"Name"], data.subset[j,"adjp"]))
        title(xlab="")
        title(ylab="Normalised Counts")
        legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
               pt.bg=1:length(levels(samples$condition)))
    }
    graphics.off()
}
