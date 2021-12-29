library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(cowplot)


#--------------all repeats-----------------------------------------------------

repeats = read.csv('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/genome_annotation_all_repeats/hvulg_rptcrft_merged.rmerge_renamed_sorted.gff',
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
  labs(title = 'Distribution of different repeat classes in Hydra vulgaris',
       x = '', y = 'Number of repeats')

ggsave(plot = all_genes_plot, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/freq_rep_all_repeats.tiff',
       width = 10, height = 10, device = 'tiff')


#--------------------------UP DEG--------------------------------------------------------

up_repeats = read.csv('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/upregulated_gene_annotation.tsv',
                      stringsAsFactors = F, header = F, sep = '\t')

up_repeats$Repeats = sapply(up_repeats$V3, function(x) {strsplit(x, '_')[[1]][1]})

freq_rep_up = as.data.frame(table(up_repeats$Repeats))

str(freq_rep_up)
freq_rep_up$Var1 = as.character(freq_rep_up$Var1)
freq_rep_up$Var1[freq_rep_up$Var1 == 'Low'] = 'Low_complexity'
freq_rep_up$Var1[freq_rep_up$Var1 == 'Simple'] = 'Simple_repeat'

up_freq_merged = merge(freq_rep_up, freq_rep_all, by = 'Var1')

up_freq_merged$Repeat_norm = up_freq_merged$Freq.x/up_freq_merged$Freq.y/length(unique(up_repeats$V19))
up_freq_merged$Repeat_class = sapply(up_freq_merged$Var1, function(x) {strsplit(x, '/')[[1]][1]})


up_plot = ggplot(up_freq_merged, aes(reorder(Var1, -Repeat_norm), Repeat_norm*100000)) +
  geom_bar(stat = 'identity', position = 'dodge', fill = 'gray88', color = 'gray66') +
  theme_bw() +
  theme(
    #plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size =18,color = "black")) +
#  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,3,4,8,9,10,11)]) +
 # scale_fill_manual(name = "Repeats", values=c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[3], "grey55")) +
  labs(
    #title = 'Distribution of repeat classes nearby upregulated genes during regeneration',
       x = 'Repeat families', y = 'Percentage of repeats') +
  scale_y_continuous(breaks=seq(0, 8, 2), limits=c(0, 8))

#---------------DOWN DEG--------------------------------------------------------------

down_repeats = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/downregulated_gene_annotation.tsv',
                      stringsAsFactors = F, header = F, sep = '\t')

down_repeats$Repeats = sapply(down_repeats$V3, function(x) {strsplit(x, '_')[[1]][1]})

freq_rep_down = as.data.frame(table(down_repeats$Repeats))

str(freq_rep_down)
freq_rep_down$Var1 = as.character(freq_rep_down$Var1)
freq_rep_down$Var1[freq_rep_down$Var1 == 'Low'] = 'Low_complexity'
freq_rep_down$Var1[freq_rep_down$Var1 == 'Simple'] = 'Simple_repeat'

down_freq_merged = merge(freq_rep_down, freq_rep_all, by = 'Var1')
down_freq_merged$Repeat_norm = down_freq_merged$Freq.x/down_freq_merged$Freq.y/length(unique(down_repeats$V19))

down_plot = ggplot(down_freq_merged, aes(reorder(Var1, -Repeat_norm), Repeat_norm*100000, fill=factor(ifelse(grepl('LINE/CR1', Var1), 'LINE/CR1', ifelse(grepl('LINE/L2', Var1), 'LINE/L2', 'Other'))))) +
  geom_bar(stat = 'identity', position = 'dodge', fill = 'gray88', color = 'gray66') +
  theme_bw() +
  theme(
    #plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size =18,color = "black")) +
#  scale_fill_manual(name = "Repeats", values=c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[3], "grey55")) +
  labs(
    #title = 'Distribution of repeat classes nearby downregulated during regeneration genes',
       x = 'Repeat families', y = 'Percentage of repeats') +
  scale_y_continuous(breaks=seq(0, 8, 2), limits=c(0, 8))

#-------------------all genes---------------------------

all_genes = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/all_repeats_gene_annotation_filtered_dist.tsv',
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
    plot.title = element_blank(),
    axis.title.x = element_text(size=18,color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=18,color = "black"),
    axis.title.y = element_text(size=18,color = "black"),
    legend.text = element_text(size=18,color = "black"),
    legend.title = element_text(size =18,color = "black")
  ) +
  labs(x = 'Repeat families', y = 'Percentage of repeats') +
  scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 104))

ggsave(plot = all_plot, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/all_genes_repeat_distribution.tiff',
       width = 10, height = 10, device = 'tiff')

#------combined plot---------
title <- ggdraw() + 
  draw_label(
    "Distribution of repeat classes nearby regeneration genes",
    fontface = 'bold', size = 18,
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(-10, 0, 0, 100)
  )

plot_row <- plot_grid(up_plot, down_plot,all_plot, labels = c('Upregulated', 'Downregulated', '        All'),
                      label_size = 18,
                      hjust = -0.5, vjust = -0.5,
                      rel_widths = c(1, 1), rel_heights = c(1,1))

combined_plots = plot_grid(title, plot_row,
                           ncol = 1,
                           rel_heights = c(0.1, 1))

ggsave(plot = combined_plots,
       filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/up_down_regeneration_distribution.tiff',
       width = 15, height = 10, device = 'tiff')


