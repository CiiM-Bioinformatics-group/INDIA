load('output/data.RData')

peak_anno$group <- NA
peak_anno[grepl(x = peak_anno$Detailed.Annotation, pattern = 'exon', fixed = T), ]$group <- 'Exon'
peak_anno[grepl(x = peak_anno$Detailed.Annotation, pattern = 'TTS', fixed = T), ]$group <- 'TTS'
peak_anno[grepl(x = peak_anno$Detailed.Annotation, pattern = 'promoter-TSS', fixed = T), ]$group <- 'promoter-TSS'
peak_anno[grepl(x = peak_anno$Detailed.Annotation, pattern = 'Intergenic', fixed = T), ]$group <- 'Intergenic'
peak_anno[grepl(x = peak_anno$Detailed.Annotation, pattern = 'intron', fixed = T), ]$group <- 'Intron'

table(peak_anno$group) %>% as.data.frame() %>% set_colnames(c('group', 'value')) -> df

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos), 
         nudge = ifelse(group %in% c('Intergenic', 'Intron'), 0, 1))

pdf('output/peak_anno.pdf', width = 5, height = 5)
ggplot(df, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = group),
                   size = 4.5, nudge_x = df2$nudge, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void() +
  theme(legend.position = 'none', plot.title = element_text(hjust = .5)) +
  labs(title  = 'Peak annotation')
dev.off()
