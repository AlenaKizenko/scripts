library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(cowplot)


#--------------all repeats-----------------------------------------------------

repeats = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_rptcrft_merged.rmerge_renamed_sorted.gff',
                   stringsAsFactors = F, header = F, sep = '\t')

repeats$Repeats = sapply(repeats$V3, function(x) {strsplit(x, '_')[[1]][1]})

freq_rep_all = as.data.frame(table(repeats$Repeats))

str(freq_rep_all)
freq_rep_all$Var1 = as.character(freq_rep_all$Var1)
freq_rep_all$Var1[freq_rep_all$Var1 == 'Low'] = 'Low_complexity'
freq_rep_all$Var1[freq_rep_all$Var1 == 'Simple'] = 'Simple_repeat'


freq_rep_all = subset(freq_rep_all, Freq > 200)
freq_rep_all = subset(freq_rep_all, Var1 != 'Unknown')
freq_rep_all = subset(freq_rep_all, Var1 != 'Low_complexity')
freq_rep_all = subset(freq_rep_all, Var1 != 'Simple_repeat')

all_genes_plot = ggplot(freq_rep_all, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = 'identity', position = 'dodge', fill = 'gray88', color = 'gray66') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=8,color = "black", angle = 40, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 8,colour = "black"),
        legend.text = element_text(size=12,color = "black"),
        legend.title = element_text(size =12,color = "black")) +
  labs(title = 'Distribution of different repeat classes in Hydra viridissima genome',
       x = '', y = 'Number of repeats')

ggsave(plot = all_genes_plot, filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/freq_rep_all_repeats_hvir.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------------all genes---------------------------

all_genes = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_repeats_near_genes_filtered.tsv',
                     stringsAsFactors = F, header = F, sep = '\t')

View(head(all_genes))
all_genes$Repeats = sapply(all_genes$V3, function(x) {strsplit(x, '_')[[1]][1]})

freq_rep = as.data.frame(table(all_genes$Repeats))

str(freq_rep)
freq_rep$Var1 = as.character(freq_rep$Var1)
freq_rep$Var1[freq_rep$Var1 == 'Low'] = 'Low_complexity'
freq_rep$Var1[freq_rep$Var1 == 'Simple'] = 'Simple_repeat'
freq_merged = merge(freq_rep, freq_rep_all, by = 'Var1')
freq_merged$Repeat_norm = freq_merged$Freq.x/freq_merged$Freq.y


all_plot = ggplot(freq_merged, aes(reorder(Var1, -Repeat_norm), Repeat_norm*100)) +
  geom_bar(stat = 'identity', position = 'dodge', fill = 'gray88', color = 'gray66') +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
    axis.title.x = element_text(size=18,color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=18,color = "black"),
    axis.title.y = element_text(size=18,color = "black"),
    legend.text = element_text(size=18,color = "black"),
    legend.title = element_text(size =18,color = "black")
  ) +
  labs(x = 'Repeat families', y = 'Percentage of repeats') +
  scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 135)) +
  labs(title = 'Percentage of repeats near genes in Hydra viridissima genome')

ggsave(plot = all_plot, filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/all_genes_repeat_distribution_hvir.tiff',
       width = 10, height = 10, device = 'tiff')

#

df_filtered = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_repeats_near_genes_filtered.tsv',
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
    x = 'Repeat families', y = 'Distance, kb')
 # scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,3,4,5,6,7,8,9,10,11)]) +
  + scale_y_continuous(breaks=seq(0, 8, 2), limits=c(0, 10))

ggsave(plot = plot_box_filt, filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/dist_repeat_distribution_hvir.tiff',
       width = 10, height = 10, device = 'tiff')
