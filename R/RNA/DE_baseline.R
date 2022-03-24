# Compare T0 RPMI Europe vs T0 RPMI India
# Baseline population effect
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

pval.thresh = 0.05
logfc.thresh = 0.0

meta.sub <- pheno %>% filter(stimulation == 'RPMI') %>% filter(time == 'Before')
meta.sub$sampleID <- factor(meta.sub$sampleID)
meta.sub$time <- factor(meta.sub$time)

counts.sub <- counts[, meta.sub$sample]

# Filter lowly expressed garbage
cutoff <- 20
table(rowSums(counts.sub) > cutoff)
counts.sub <- counts.sub[rowSums(counts.sub) > cutoff, ]

dds <- DESeqDataSetFromMatrix(countData = counts.sub,
                            colData = meta.sub,
                            design = ~ population + sex) # TODO: Add gender / age when available

dds <- DESeq(dds)
res <- results(dds, contrast = c('population', 'Europe', 'India')) %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

# merge with the amount of reads per group. Median
samples1 <- meta.sub %>% filter(population == 'Europe') %>% pull(sample)
median1 = counts.sub[, samples1] %>% as.matrix() %>% rowMedians()

samples2 <- meta.sub %>% filter(population == 'India') %>% pull(sample)
median2 = counts.sub[, samples2] %>% as.matrix() %>% rowMedians()

res <- cbind(res,
           data.frame(
             'median_reads_europe' = median1,
             'median_reads_india' = median2
           )) %>% 
arrange(padj) %>% 
mutate(significance = ifelse(padj < pval.thresh & abs(log2FoldChange) > logfc.thresh, TRUE, FALSE),
       direction = ifelse(log2FoldChange > 0, 'Upregulated', 'Downregulated'))

write.xlsx(x = res, 'output/DE/BL/europe_T0RPMI_vs_india_T0RPMI.xlsx', overwrite=T)







