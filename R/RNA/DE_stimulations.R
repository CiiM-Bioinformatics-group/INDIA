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

runDEseq2 <- function(counts, meta, pop, control, stimul, t, outfile, pval.thresh, logfc.thresh) {

  meta.sub <- meta %>%
    filter(population == pop) %>%
    filter(stimulation == stimul | stimulation == control) %>%
    filter(time == t)

  meta.sub$sampleID <- factor(meta.sub$sampleID)
  print(table(meta.sub$stimulation, meta.sub$time))

  counts.sub <- counts[, meta.sub$sample]

  # Filter lowly expressed garbage
  cutoff <- 20
  table(rowSums(counts.sub) > cutoff)
  counts.sub <- counts.sub[rowSums(counts.sub) > cutoff, ]

  dds <- DESeqDataSetFromMatrix(countData = counts.sub,
                                colData = meta.sub,
                                design = ~ stimulation + sampleID) # TODO: Add gender / age when available

  dds <- DESeq(dds)

  res <- results(dds, contrast = c('stimulation', stimul, control)) %>%
    as.data.frame() %>%
    # arrange(padj) %>%
    mutate(gene = rownames(.))

  # merge with the amount of reads per group. Median
  samples1 <- meta.sub %>% filter(stimulation == stimul) %>% pull(sample)
  median1 = counts.sub[, samples1] %>% as.matrix() %>% rowMedians()

  samples2 <- meta.sub %>% filter(stimulation == control) %>% pull(sample)
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
  # India
  c(pop = 'India', control = 'RPMI', stimul = 'SARS-CoV-2', t = 'Before', outfile = 'output/DE/stimulations/india_T0RPMI_vs_T0SARS.xlsx'),
  c(pop = 'India', control = 'RPMI', stimul = 'Influenza', t = 'Before', outfile = 'output/DE/stimulations/india_T0RPMI_vs_T0INFL.xlsx'),
  c(pop = 'India', control = 'RPMI', stimul = 'SARS-CoV-2', t = 'After', outfile = 'output/DE/stimulations/india_T2RPMI_vs_T2SARS.xlsx'),
  c(pop = 'India', control = 'RPMI', stimul = 'Influenza', t = 'After', outfile = 'output/DE/stimulations/india_T2RPMI_vs_T2INFL.xlsx'),

  # Europe
  c(pop = 'Europe', control = 'RPMI', stimul = 'SARS-CoV-2', t = 'Before', outfile = 'output/DE/stimulations/europe_T0RPMI_vs_T0SARS.xlsx'),
  c(pop = 'Europe', control = 'RPMI', stimul = 'Influenza', t = 'Before', outfile = 'output/DE/stimulations/europe_T0RPMI_vs_T0INFL.xlsx'),
  c(pop = 'Europe', control = 'RPMI', stimul = 'SARS-CoV-2', t = 'After', outfile = 'output/DE/stimulations/europe_T2RPMI_vs_T2SARS.xlsx'),
  c(pop = 'Europe', control = 'RPMI', stimul = 'Influenza', t = 'After', outfile = 'output/DE/stimulations/europe_T2RPMI_vs_T2INFL.xlsx')
)

for (comp in comps) {
  runDEseq2(counts = counts,
            meta = pheno,
            pop = comp[['pop']],
            control = comp[['control']],
            stimul = comp[['stimul']],
            t = comp[['t']],
            pval.thresh = PVAL.THRESH,
            logfc.thresh = LOGFC.THRESH,
            outfile = comp[['outfile']]
          )
}
