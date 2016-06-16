#!/usr/bin/env Rscript

library(derfinder)

Args      <- commandArgs()
chr       <- ifelse(is.na(Args[6]),  NA,            Args[6])
cutoff    <- ifelse(is.na(Args[7]),  10, as.integer(Args[7]))
readLen   <- ifelse(is.na(Args[8]),  75, as.integer(Args[8]))
bigWigDir <- ifelse(is.na(Args[9]),  "tophat",      Args[9])
chrFile   <- ifelse(is.na(Args[10]), "chrom.sizes", Args[10])
species   <- ifelse(is.na(Args[11]), "danio_rerio", Args[11])

options(species=species)
options(chrsStyle='NCBI')

files <- rawFiles(datadir=bigWigDir, samplepatt='*.bw$', fileterm=NULL)
names(files) <- gsub('.bw', '', names(files))

chr.sizes <- read.table(chrFile)
chr.sizes$V1 <- as.character(chr.sizes$V1)
chrs <- chr.sizes$V1[chr.sizes$V1 == chr]
chrlens <- chr.sizes$V2[chr.sizes$V1 == chr]
if (chr == "spike") {
    chrs <- chr.sizes$V1[grepl("ERCC", chr.sizes$V1)]
    chrlens <- chr.sizes$V21[grepl("ERCC", chr.sizes$V1)]
}

fullCov <- fullCoverage(files=files, chrs=chrs, chrlens=chrlens, cutoff=cutoff)

regionMat <- regionMatrix(fullCov, cutoff=cutoff, L=readLen)

write.table(regionMat[[1]]$regions, file=paste(chr, ".regions", sep=""),
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(regionMat[[1]]$coverageMatrix, file=paste(chr, ".matrix", sep=""),
            quote=FALSE, sep="\t", row.names=FALSE)
