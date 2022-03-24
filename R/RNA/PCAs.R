try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(DESeq2)

setwd('/vol/projects/CIIM/INDIA/analysis/RNA')
load('output/data.Rdata')

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = pheno, 
                              design = ~ population + time + stimulation + sex)
vsd <- varianceStabilizingTransformation(dds)
pca <- plotPCA(vsd, intgroup = c('population', 'time', 'stimulation', 'sex'), returnData=T)

pdf('output/pcas.pdf', width = 16, height = 4)
cowplot::plot_grid(
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = stimulation)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_nejm() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Stimulation'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = population)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_uchicago() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Population'),
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = time)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_jama() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Time'), 
  
  ggplot(pca) +
    geom_point(aes(PC1, PC2, color = sex)) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_viridis_d() +
    labs(x = 'PC1 [47%]', y = 'PC2 [17%]', color = ' ', title = 'Sex'), 
  
  
  ncol = 4, nrow = 1, align = 'hv', axis = 'tblr'
  
)
dev.off()




# Per population
for (p in c('India', 'Europe')) {
  
  meta.sub <- pheno %>% filter(population == p)
  counts.sub <- counts[, meta.sub$sample]
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts.sub, 
                                colData = meta.sub, 
                                design = ~ time + stimulation + sex)
  vsd <- varianceStabilizingTransformation(dds)
  
  tmp <- plotPCA(vsd, intgroup = c('time'))
  xlab <- tmp$labels$x
  ylab <- tmp$labels$y
  
  pca <- plotPCA(vsd, intgroup = c('time', 'stimulation', 'sex'), returnData=T)
  
  pdf(paste0('output/PCA_pop', p, '.pdf'), width = 12, height = 4)
  print(cowplot::plot_grid(
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = stimulation)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_nejm() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Stimulation'),
    
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = time)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_jama() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Time'), 
    
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = sex)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_viridis_d() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Sex'), 
    
    ncol = 3, nrow = 1, align = 'hv', axis = 'tblr'
    
  )
  )
  dev.off()
  
}


# Per stimulation
for (s in c('RPMI', 'Influenza', 'SARS-CoV-2')) {
  
  meta.sub <- pheno %>% filter(stimulation == s)
  counts.sub <- counts[, meta.sub$sample]
  
  dds <- DESeqDataSetFromMatrix(countData = counts.sub, 
                                colData = meta.sub, 
                                design = ~ population + time + sex)
  vsd <- varianceStabilizingTransformation(dds)
  
  tmp <- plotPCA(vsd, intgroup = c('time'))
  xlab <- tmp$labels$x
  ylab <- tmp$labels$y
  
  pca <- plotPCA(vsd, intgroup = c('time', 'population', 'sex'), returnData=T)
  
  pdf(paste0('output/PCA_stim_', s, '.pdf'), width = 8, height = 3)
  print(cowplot::plot_grid(
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = population)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_uchicago() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Population'),
    
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = time)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_jama() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Time'), 
    
    ggplot(pca) +
      geom_point(aes(PC1, PC2, color = sex)) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      scale_color_viridis_d() +
      labs(x = xlab, y = ylab, color = ' ', title = 'Sex'),
    
    
    ncol = 3, nrow = 1, align = 'hv', axis = 'tblr'
    
    )
  )
  dev.off()
  
}
