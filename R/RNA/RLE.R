# RLE plot to check for outliers in bulkRNAseq

library(ggplot2)
library(dplyr)
library(magrittr)
library(reshape2)
library(DESeq2)
library(ggsci)
library(EDASeq)

rm(list = ls())
setwd('/vol/projects/CIIM/INDIA/analysis/RNA/')
load('output/data.Rdata')

# Re-order data in the same way as multiQC for consistency
df <- read.csv('results/multiqc/star/multiqc_data/mqc_picard_deduplication_1.txt', sep = '\t', row.names = 1)
counts <- counts[, rownames(df)]

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = pheno, 
                              design = ~ population + time + stimulation)
dds <- DESeq(dds)

plotRLE(DESeq2::counts(dds, normalized = TRUE)) -> x

x %>% 
  as.data.frame() %>% 
  mutate(gene = rownames(.)) %>% 
  reshape2::melt(., id.vars = 'gene') -> 
  y

# Figure out the samples that need to be colored
samples <- c("india_T2_SARS_R4", 'india_T2_RPMI_R4','india_T0_INFL_R6', 'europe_T2_SARS_R6')

y$color <- NA
y[which(y$variable %in% samples), ]$color <- 'Warning'

pdf('output/RLE.pdf', width = 5, height = 6)
ggplot(y) +
  geom_boxplot(aes(x = value, y = variable, fill = color), outlier.shape = NA) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_classic() +
  labs(x = 'RLE', y = 'Samples', title = 'Relative log expression') +
  theme(plot.title = element_text(hjust = .5), axis.text.y = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values = c('red'))
dev.off()
