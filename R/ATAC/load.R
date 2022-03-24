try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/ATAC')

# Count matrix
counts <- read.csv('results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt', sep = '\t', comment.char = '#')
rownames(counts) <- counts$Geneid
counts %<>% select(-all_of(c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length')))
colnames(counts) <- lapply(strsplit(x = colnames(counts), split = '.', fixed = T), function(x){x[[1]]} ) 


# Peak annotation
peak_anno <- read.csv('results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt', sep = '\t', row.names = 1)
peak_anno <- peak_anno[rownames(counts), ]

# Phenotypes
pheno <- read.csv('results/pipeline_info/design_reads.csv', header=T, sep = ',')
pheno$file <- lapply(strsplit(x = pheno$fastq_1, split = '/'), function(x){x[[8]]} ) %>% unlist()
pheno$sample <- substr(pheno$sample, start = 0, stop = nchar(pheno$sample) - 3)

pheno <- cbind(pheno, 
               colsplit(pheno$file, pattern = '_', names = c('population', 'sampleID', 'time', 'trash'))) %>% 
  select(-trash) %>% 
  select(all_of(c('sample', 'population', 'sampleID', 'time')))

pheno$time <- ifelse(pheno$time == 'T0', 'Before', 'After')

pheno %<>% arrange(match(sample, colnames(counts)))

# Add genders and ages
meta <- read.csv('../metadata.csv', sep = ';', header=T)
meta$sex <- ifelse(meta$sex == 'Male', 'male', 'female')
meta$population <- ifelse(meta$population == 'Europe', 'europe', 'india')

pheno <- merge(pheno, meta, by = c('population', 'sampleID')) %>% arrange(match(sample, colnames(counts)))


# QC and filtering
samples.rem <- c('europe_T2_R7', 'india_T2_R1')
peaks.rem <- c()

pheno %<>% filter(!sample %in% samples.rem)
counts <- counts[, pheno$sample]

stopifnot(all(pheno$sample == colnames(counts)))
stopifnot(all(rownames(counts) == rownames(counts)))

save.image('output/data.RData')





