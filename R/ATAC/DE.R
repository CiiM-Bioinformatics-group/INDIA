# Differences in chromatin accessibility between European and Indian people
# at the baseline level.

try(dev.off())
rm(list = ls())
library(DESeq2)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(openxlsx)
setwd('/vol/projects/CIIM/INDIA/analysis/ATAC')
load('output/data.RData')

pval.thresh = 0.05
logfc.thresh = 0.0


# Remove the X and Y peaks -> Not interested in sex effect
head(peak_anno)
peak_anno %<>% filter(Chr != 'chrX') %>% filter(Chr != "chrY")
counts <- counts[rownames(peak_anno), ]

runDEseq2 <- function(counts, meta, t, outfile, pval.thresh, logfc.thresh) {
  
  meta.sub <- meta %>%
    filter(time == t)
  
  meta.sub$sampleID <- factor(meta.sub$sampleID)
  print(meta.sub)
  counts.sub <- counts[, meta.sub$sample]
  
  # Filter lowly expressed garbage
  cutoff <- 20
  table(rowSums(counts.sub) > cutoff)
  counts.sub <- counts.sub[rowSums(counts.sub) > cutoff, ]
  
  dds <- DESeqDataSetFromMatrix(countData = counts.sub,
                                colData = meta.sub,
                                design = ~ population + sex)
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c('population', 'europe', 'india')) %>%
    as.data.frame() %>%
    # arrange(padj) %>%
    mutate(gene = rownames(.))
  
  # merge with the amount of reads per group. Median
  samples1 <- meta.sub %>% filter(population == 'europe') %>% pull(sample)
  median1 = counts.sub[, samples1] %>% as.matrix() %>% rowMedians()
  
  samples2 <- meta.sub %>% filter(population == 'india') %>% pull(sample)
  median2 = counts.sub[, samples2] %>% as.matrix() %>% rowMedians()
  
  res <- cbind(res,
               data.frame(
                 'median_reads_stimul' = median1,
                 'median_reads_control' = median2
               )) %>%
    arrange(padj) %>%
    mutate(significance = ifelse(padj < pval.thresh & abs(log2FoldChange) > logfc.thresh, TRUE, FALSE),
           direction = ifelse(log2FoldChange > 0, 'Upregulated', 'Downregulated'))
  
  write.xlsx(res, outfile, overwrite = T)
}


comps <- list(
  c(t = 'Before', outfile = 'output/DE/europe_vs_india_before.xlsx'),
  c(t = 'After', outfile = 'output/DE/europe_vs_india_after.xlsx')
)


for (comp in comps) {
  runDEseq2(counts = counts,
            meta = pheno,
            t = comp[['t']],
            pval.thresh = pval.thresh,
            logfc.thresh = logfc.thresh,
            outfile = comp[['outfile']]
  )
}  
