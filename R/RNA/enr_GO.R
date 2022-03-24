try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(ggsci)

setwd('/vol/projects/CIIM/INDIA/analysis/RNA')
n.show <- 5
n.thresh <- 5

deg <- list(
  # India
  'India SARS-CoV-2 Before' = read.xlsx('output/DE/stimulations/india_T0RPMI_vs_T0SARS.xlsx'), 
  'India Influenza Before' = read.xlsx('output/DE/stimulations/india_T0RPMI_vs_T0INFL.xlsx'),
  
  'India SARS-CoV-2 After' = read.xlsx('output/DE/stimulations/india_T2RPMI_vs_T2SARS.xlsx'),
  'India Influenza After' = read.xlsx('output/DE/stimulations/india_T2RPMI_vs_T2INFL.xlsx'),
  
  # Europe
  'Europe SARS-CoV-2 Before' = read.xlsx('output/DE/stimulations/europe_T0RPMI_vs_T0SARS.xlsx'), 
  'Europe Influenza Before' = read.xlsx('output/DE/stimulations/europe_T0RPMI_vs_T0INFL.xlsx'),
  
  'Europe SARS-CoV-2 After' = read.xlsx('output/DE/stimulations/europe_T2RPMI_vs_T2SARS.xlsx'),
  'Europe Influenza After' = read.xlsx('output/DE/stimulations/europe_T2RPMI_vs_T2INFL.xlsx')
)


enrich <- function(genes, thresh.incl=5, n.show=5) {
  
  if (length(genes) == 0) {
    print('No up/downregulated genes. Skipping enrichment')
    return()
  }
  
  enr.GO <- enrichGO(gene = genes,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  return(enr.GO)
}
getUp <- function(x) {x %>% filter(significance == T) %>% filter(direction == 'Upregulated') %>% pull(gene)}
getDown <- function(x) {x %>% filter(significance == T) %>% filter(direction == 'Downregulated') %>% pull(gene)}
conv <- function(x) { bitr(geneID = x, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb="org.Hs.eg.db") %>% pull(ENTREZID) }

# Get genes, convert to Entrez IDs, run enrichment and merge result for down and upregulated separately
deg %>% lapply(., getUp) %>% lapply(., conv) %>% lapply(., enrich) %>% merge_result() -> up
deg %>% lapply(., getDown) %>% lapply(., conv) %>% lapply(., enrich) %>% merge_result() -> down

write.xlsx(x = up, file = 'output/enr_res_up_stimulations.xlsx')
write.xlsx(x = down, file = 'output/enr_res_down_stimulations.xlsx')

# Upreg.
df <- up
df <- data.frame(df@compareClusterResult)

df <- cbind(df, 
            colsplit(df$Cluster, pattern = ' ', names = c('Population', 'Stimulation', 'Time'))) %>% 
  mutate(GeneRatio2 = sapply(df$GeneRatio, function(x) eval(parse(text=x)))) %>%
  filter(Count >= n.thresh) %>% 
  # mutate(Description = stringr::str_wrap(.$Description, width = 50)) %>%
  group_by(Cluster) %>%
  dplyr::slice(1:n.show)

# Order everything in the way that we want
df$Time <- factor(df$Time, levels = c('Before', 'After'))
df$Stimulation <- factor(df$Stimulation, levels = c('SARS-CoV-2', 'Influenza'))
df$Description <- factor(df$Description, levels = unique(df$Description))

pdf('output/enr_up.pdf', width = 12, height = 10)
print(ggplot(data = df) +
        geom_point(aes(x = Time, y = Description, size = GeneRatio2, color = p.adjust)) +
        DOSE::theme_dose() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
              strip.background = element_rect(fill = 'white', color = 'black'),
              strip.text = element_text(face = 'bold', size = 15),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(), 
              plot.title = element_text(hjust = 0.5), 
        aspect.ratio = nrow(df) * 0.075) +
        scale_colour_gradient(limits=c(0, 0.05), low="darkred", high = 'darkblue') +
        scale_x_discrete(drop=F) +
        labs(color = 'Adj. P', size = 'Gene ratio', title = 'Upregulated processes') +
        facet_grid(Stimulation ~ Population)
      )
dev.off()




# Downreg
df <- down
df <- data.frame(df@compareClusterResult)

df <- cbind(df, 
            colsplit(df$Cluster, pattern = ' ', names = c('Population', 'Stimulation', 'Time'))) %>% 
  mutate(GeneRatio2 = sapply(df$GeneRatio, function(x) eval(parse(text=x)))) %>%
  filter(Count >= n.thresh) %>% 
  # mutate(Description = stringr::str_wrap(.$Description, width = 50)) %>%
  group_by(Cluster) %>%
  dplyr::slice(1:n.show)

# Order everything in the way that we want
df$Time <- factor(df$Time, levels = c('Before', 'After'))
df$Stimulation <- factor(df$Stimulation, levels = c('SARS-CoV-2', 'Influenza'))
df$Description <- factor(df$Description, levels = unique(df$Description))

pdf('output/enr_down.pdf', width = 12, height = 10)
print(ggplot(data = df) +
        geom_point(aes(x = Time, y = Description, size = GeneRatio2, color = p.adjust)) +
        DOSE::theme_dose() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
              strip.background = element_rect(fill = 'white', color = 'black'),
              strip.text = element_text(face = 'bold', size = 15),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              aspect.ratio = nrow(df) * 0.075) +
        scale_colour_gradient(limits=c(0, 0.05), low="darkred", high = 'darkblue') +
        scale_x_discrete(drop=F) +
        labs(color = 'Adj. P', size = 'Gene ratio', title = 'Downregulated processes') +
        facet_grid(Stimulation ~ Population)
)
dev.off()
