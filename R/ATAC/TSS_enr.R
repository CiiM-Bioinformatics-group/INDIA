rm(list = ls())
try(dev.off())

library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(parallel)

setwd('/vol/projects/CIIM/INDIA/analysis/ATAC/')

dir <- "/vol/projects/CIIM/INDIA/analysis/ATAC/results/bwa/mergedLibrary/"
bam.files <- list.files(dir, pattern='*bam$')

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

calc <- function(bamfile, txs) {
  bam <- readBamFile(bamFile = paste0(dir, bamfile),
                     bigFile = T)

  x <- TSSEscore(obj = bam, txs = txs)
  
  print(paste0(bamfile, ':', x))
  res <- list()
  res[[bamfile]] <- x$values

  return (res)

}

res <- mclapply(X = bam.files, FUN=calc, txs = txs, mc.cores=10)
print(res)
df <- do.call(cbind.data.frame, res)
write.csv(df, file = 'output/TSSenrichment.csv')