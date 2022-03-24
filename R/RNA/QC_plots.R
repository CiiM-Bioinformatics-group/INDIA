# QC plots for the bulkRNAseq: duplication levels per sample and STAR mapping overview.

try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/RNA')

# PICARD markduplicates
df <- read.csv('results/multiqc/star/multiqc_data/mqc_picard_deduplication_1.txt', sep = '\t', row.names = 1) %>% 
  set_colnames(c('Uniquely pairs', 'Unique unpaired', 'Duplicate pairs (optical)', 'Duplicate pairs (non-optical)', 'Duplicate unpaired')) 

# Calc the percentages
sums <- apply(X = df, MARGIN = 1, FUN = sum)
df <- apply(X = df, MARGIN = 2, FUN = function(x, sums) {round(x / sums, 2) * 100}, sums = sums) %>% 
  as.data.frame() %>% 
  mutate(sampleID = rownames(.)) %>% 
  melt()

dups <- ggplot(df) +
  geom_bar(aes(x = value, y = sampleID, fill = variable), position = 'stack', stat = 'identity') +
  theme_classic() +
  scale_fill_npg() +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = '% of reads', y = 'Samples', title = 'Duplication levels per sample', fill = 'Category')

dups

# STAR mapping
df <- read.csv('results/multiqc/star/multiqc_data/mqc_star_alignment_plot_1.txt', sep = '\t', row.names = 1) %>% 
  mutate(sampleID = rownames(.)) %>% 
  set_colnames(c('Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped (too short)', 'Unmapped (other)', 'sampleID')) %>% 
  melt()

df$variable <- factor(df$variable, 
                       levels = c('Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped (too short)', 'Unmapped (other)'))

nreads <- ggplot(df) +
  geom_bar(aes(x = value, y = sampleID, fill = variable), position = 'stack', stat = 'identity') +
  theme_classic() +
  scale_fill_npg() +
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Nr. reads', y = 'Samples', title = 'Total nr reads per sample', fill = 'Category')


# Combine and output
pdf('output/qc_plots.pdf', width = 10, height = 6)
ggpubr::ggarrange(dups, nreads)
dev.off()

# median nr of reads
df <- read.csv('results/multiqc/star/multiqc_data/multiqc_featureCounts.txt', sep = '\t', row.names=1)
median(df$Total) / 1e6

