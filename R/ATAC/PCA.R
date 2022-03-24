try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/ATAC')
load('output/data.RData')

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = pheno, 
                              design = ~ population + time + sex)

vsd <- varianceStabilizingTransformation(dds)

pca <- plotPCA(vsd, intgroup = c('population', 'time', 'sex'), returnData=T)

pdf('output/pcas.pdf', width = 10, height = 4)
ggpubr::ggarrange(
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = population)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_uchicago() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Population'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = sex)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_colour_tron() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Sex'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = time)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_jama() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Time'), ncol = 3, nrow = 1
  
)
dev.off()



# Wthout X chromosome
peak_anno %<>% filter(Chr != 'chrX') %>% filter(Chr != 'chrY')
counts <- counts[rownames(peak_anno), ]


dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = pheno, 
                              design = ~ population + time + sex)

vsd <- varianceStabilizingTransformation(dds)
pca <- plotPCA(vsd, intgroup = c('population', 'time', 'sex'), returnData=T)

pdf('output/pcas_nogender.pdf', width = 10, height = 4)
ggpubr::ggarrange(
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = population)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_uchicago() +
    labs(x = 'PC1 [13%]', y = 'PC2 [11%]', color = ' ', title = 'Population'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = sex)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_colour_tron() +
    labs(x = 'PC1 [13%]', y = 'PC2 [11%]', color = ' ', title = 'Sex'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = time)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_jama() +
    labs(x = 'PC1 [13%]', y = 'PC2 [11%]', color = ' ', title = 'Time'), ncol = 3, nrow = 1
  
)
dev.off()
