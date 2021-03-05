library(ggplot2)
library(RColorBrewer)
library(cowplot)


df = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/all_repeats_gene_annotation_dist.tsv',
              stringsAsFactors = F, header = F, sep = '\t')

View(head(df))
df$Repeat = sapply(df$V2, function(x) {strsplit(x, '_')[[1]][1]})
df$Repeat[df$Repeat == 'Low'] = 'Low complexity'
df$Repeat[df$Repeat == 'Simple'] = 'Simple repeats'
df$Repeat_class = sapply(df$Repeat, function(x) {strsplit(x, '/')[[1]][1]})
str(df)
df$Repeat_class = as.factor(df$Repeat_class)
plot_box = ggplot(df, aes(Repeat, V8/1000, fill = Repeat_class)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, alpha = 0.8) +
  theme_bw() +
  theme(
    #plot.title = element_text(hjust = 0, size = 16, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none',
       # legend.title = element_text(size = 14),
        #legend.text = element_text(size = 14),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black")
  ) +
  labs(
    #title = 'Distances distribution from various repeat families to nearest genes',
       x = 'Repeat families', y = 'Distance, kb') +
  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,3,4,8,9,10,11)])

#------after filtration by 10kb---------

df_filtered = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/all_repeats_gene_annotation_filtered_dist.tsv',
              stringsAsFactors = F, header = F, sep = '\t')

View(head(df_filtered))
df_filtered$Repeat = sapply(df_filtered$V3, function(x) {strsplit(x, '_')[[1]][1]})
df_filtered$Repeat[df_filtered$Repeat == 'Low'] = 'Low complexity'
df_filtered$Repeat[df_filtered$Repeat == 'Simple'] = 'Simple repeats'
df_filtered$Repeat_class = sapply(df_filtered$Repeat, function(x) {strsplit(x, '/')[[1]][1]})
str(df_filtered)
df_filtered$Repeat_class = as.factor(df_filtered$Repeat_class)
plot_box_filt = ggplot(df_filtered, aes(Repeat, V19/1000, fill = Repeat_class)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, alpha = 0.8) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 16, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black")
  ) +
  labs(
    #title = 'Distances distribution from various repeat families to nearest genes',
       x = 'Repeat families', y = 'Distance, kb') +
  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,3,4,8,9,10,11)]) +
  scale_y_continuous(breaks=seq(0, 8, 2), limits=c(0, 10))

title <- ggdraw() + 
  draw_label(
    "Distances distribution from various repeat families to nearest genes",
    fontface = 'bold', size = 18,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 380)
  )

plot_row <- plot_grid(plot_box, plot_box_filt, labels = c('                           Before filtering',
                                                          '                            After filtering'),
                      label_size = 18,
                      #label_x = 0, label_y = 0,
                      hjust = -0.5, vjust = -0.5,
                      rel_widths = c(0.8, 1), rel_heights = c(1,1))

combined_plots = plot_grid(title, plot_row,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1))

ggsave(plot = combined_plots,
       filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/combined_plots.tiff',
       width = 20, height = 10, device = 'tiff')
