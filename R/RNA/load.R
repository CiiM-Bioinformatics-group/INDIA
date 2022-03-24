library(reshape2)
library(dplyr)
library(magrittr)

rm(list = ls())
setwd('/vol/projects/CIIM/INDIA/analysis/RNA/')

# Count matrix
counts <- read.csv('results/star/featurecounts.merged.counts.tsv', sep = '\t', header = T, row.names = 1) %>% select(-gene_name)

# Phenotypes
pheno <- read.csv('results/pipeline_info/samplesheet.valid.csv', header=T, sep = ',')
pheno$file <- lapply(strsplit(x = pheno$fastq_1, split = '/'), function(x){x[[8]]} ) %>% unlist()
pheno$sample <- substr(pheno$sample, start = 0, stop = nchar(pheno$sample) - 3)
pheno <- cbind(pheno, 
               colsplit(pheno$file, pattern = '_', names = c('population', 'sampleID', 'time', 'stimulation', 'trash'))) %>% 
  select(-trash) %>% 
  select(all_of(c('sample', 'population', 'sampleID', 'time', 'stimulation')))

# Adjust some of the variables in pheno
pheno$time <- ifelse(pheno$time == 'T0', 'Before', 'After')
pheno$population <- ifelse(pheno$population == 'india', 'India', 'Europe')
pheno$stimulation <- ifelse(pheno$stimulation == 'INFL', 'Influenza', 
                            ifelse(pheno$stimulation == 'SARS', 'SARS-CoV-2', 'RPMI'))

# Add genders and ages
meta <- read.csv('../metadata.csv', sep = ';', header=T)
pheno <- merge(pheno, meta, by = c('population', 'sampleID'))


# Kick out the samples of low quality based on QC plots (QC_plots.R) and RLE plot (RLE.R)
rem <- c("india_T2_SARS_R4", 'india_T2_RPMI_R4','india_T0_INFL_R6', 'europe_T2_SARS_R6', 'europe_T0_SARS_R9')
pheno %<>% filter(!sample %in% rem)

# Kick out the 0 var genes
vars <- apply(X = counts, MARGIN = 1, FUN = var)
keep <- which(vars != 0)


counts <- counts[keep, pheno$sample]
stopifnot(all(pheno$sample == colnames(counts)))

rm(keep, vars, rem, meta)
save.image('output/data.Rdata')
