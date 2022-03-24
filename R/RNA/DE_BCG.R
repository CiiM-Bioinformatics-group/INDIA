# Compare T2 RPMI vs T0 RPMI for european and indian populations separately
# BCG vaccination effect

try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(openxlsx)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/RNA')
load('output/data.Rdata')

PVAL.THRESH = 0.05
LOGFC.THRESH = 0.0

runDEseq2 <- function(counts, meta, pop, outfile, pval.thresh, logfc.thresh) {
  
  meta.sub <- meta %>%
    filter(population == pop) %>%
    filter(stimulation == 'RPMI')

  meta.sub$sampleID <- factor(meta.sub$sampleID)
  meta.sub$time <- factor(meta.sub$time)

  print(table(meta.sub$stimulation, meta.sub$time, meta.sub$population))

  counts.sub <- counts[, meta.sub$sample]
  print(counts.sub[1:5, 1:5])
  # Filter lowly expressed garbage
  cutoff <- 20
  table(rowSums(counts.sub) > cutoff)
  counts.sub <- counts.sub[rowSums(counts.sub) > cutoff, ]

  dds <- DESeqDataSetFromMatrix(countData = counts.sub,
                                colData = meta.sub,
                                design = ~ time + sampleID)

  dds <- DESeq(dds)
  res <- results(dds, contrast = c('time', 'After', 'Before')) %>%
    as.data.frame() %>%
    mutate(gene = rownames(.))
  
  # merge with the amount of reads per group. Median
  t0.samples <- meta.sub %>% filter(time == 'Before') %>% pull(sample)
  median.t0 = counts.sub[, t0.samples] %>% as.matrix() %>% rowMedians()

  t2.samples <- meta.sub %>% filter(time == 'After') %>% pull(sample)
  median.t2 = counts.sub[, t2.samples] %>% as.matrix() %>% rowMedians()
  
  res <- cbind(res,
               data.frame(
                 'median_reads_T0' = median.t0,
                 'median_reads_T2' = median.t2
               )) %>% 
    arrange(padj) %>% 
    mutate(significance = ifelse(padj < pval.thresh & abs(log2FoldChange) > logfc.thresh, TRUE, FALSE),
           direction = ifelse(log2FoldChange > 0, 'Upregulated', 'Downregulated'))
  
  write.xlsx(res, outfile, overwrite = T)
}

comps <- list(
  list(pop = 'India', outfile = 'output/DE/BCG/india_T2RPMI_vs_T0RPMI.xlsx', pval.thresh = PVAL.THRESH, logfc.thresh = LOGFC.THRESH),
  list(pop = 'Europe', outfile = 'output/DE/BCG/europe_T2RPMI_vs_T0RPMI.xlsx', pval.thresh = PVAL.THRESH, logfc.thresh = LOGFC.THRESH)
)

for (comp in comps) {
  runDEseq2(counts = counts,
            meta = pheno,
            pop = comp[['pop']],
            pval.thresh = comp[['pval.thresh']],
            logfc.thresh = comp[['logfc.thresh']],
            outfile = comp[['outfile']]
          )
}
