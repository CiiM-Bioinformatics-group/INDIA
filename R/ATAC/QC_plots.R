try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/ATAC')


# TSS enr
df <- read.csv('output/TSSenrichment.csv', header=T, row.names=1)

colnames(df) <- colsplit(string = colnames(df), pattern = '\\.', names = c('sample', 'garbage')) %>% pull(sample)
df$distance <- 100*(-9:10-.5)

df <- melt(df, id.vars = 'distance')

pdf('output/TSS_enr.pdf', width = 4, height = 4)
ggplot(data = df) +
  geom_hline(yintercept = 5, lty=2, alpha=0.4) +
  geom_point(aes(x = distance, y = value), alpha = 0.4) +
  geom_line(aes(x = distance, y = value, group = variable), alpha = 0.4) +
  theme_classic() +
  labs(x = 'Distance to TSS', 
       y = 'TSS enrichment')
dev.off()


# Fragment size dis
# File manually downloaded from MultiQC report
df <- read.csv('output/picard_insert_size.tsv', header=T, row.names = 1, sep = '\t') %>% 
  mutate(size = as.numeric(rownames(.))) %>% 
  melt(., id.vars = 'size')

pdf('output/insert_sizes.pdf', width = 4, height = 4)
ggplot(df) +
  geom_line(aes(x = size, y = value, group = variable)) +
  theme_classic() +
  theme(legend.position = 'none', plot.title = element_text(hjust = .5)) +
  labs(x = 'Size', y = 'Nr. reads', title = 'Fragment sizes distribution per sample') +
  xlim(c(0,750))
dev.off()

# PICARD Total reads
df <- read.csv('output/picard_deduplication.tsv', sep = '\t', row.names = 1) %>% 
  set_colnames(c('Unique pairs', 'Unique unpaired', 'Duplicate pairs (optical)', 'Duplicate pairs (non-optical)', 'Duplicate unpaired', 'Unmapped')) %>% 
  # apply(., MARGIN = 1, sum) %>% as.data.frame(count = .) %>% 
  mutate(sampleID = rownames(.)) %>% 
  melt()

nreads <- ggplot(df) +
  geom_bar(aes(x = value, y = sampleID), position = 'stack', stat = 'identity') +
  theme_classic() +
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Nr. reads', y = 'Samples', title = 'Nr. reads per sample')


# PICARD markduplicates
# Manu. downloaded from multiqc report
df <- read.csv('output/picard_deduplication.tsv', sep = '\t', row.names = 1) %>% 
  set_colnames(c('Unique pairs', 'Unique unpaired', 'Duplicate pairs (optical)', 'Duplicate pairs (non-optical)', 'Duplicate unpaired', 'Unmapped'))
  
# Calc percentages
sums <- apply(X = df, MARGIN = 1, FUN = sum)
df <- apply(X = df, MARGIN = 2, FUN = function(x, sums) {round(x / sums, 5) * 100}, sums = sums) %>% 
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



# Annotation
df <- read.csv('output/map_results_multiqc.tsv', sep = '\t', row.names = 1) %>% 
  set_colnames(c("Exon", 'Intergenic', 'Intron', 'Promotor TSS', 'TTS', 'Unassigned')) %>% 
  as.data.frame() %>% 
  mutate(sampleID = rownames(.)) %>% 
  melt()

map <- ggplot(df) +
  geom_bar(aes(x = value, y = sampleID, fill = variable), position = 'stack', stat = 'identity') +
  theme_classic() +
  scale_fill_npg() +
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Nr. peaks', y = 'Samples', title = 'Total peaks per sample', fill = 'Category')

# Combine and output
pdf('output/qc_plots.pdf', width = 12, height = 6)
ggpubr::ggarrange(nreads, dups, map, nrow = 1)
dev.off()

