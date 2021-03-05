library(dplyr)
library(ggplot2)
library(topGO)
library(RColorBrewer)
library(preprocessCore)
library(ggfortify)
library(pheatmap)
library(ggsci)
library(tidyr)

#------------reading parsed data with DEG---------------------
up = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/up_parsed.tsv',
              header = F, stringsAsFactors = F, sep = '\t')

colnames(up) = c('Gene', 'Time', 'DE')

down = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/down_parsed.tsv',
                header = F, stringsAsFactors = F, sep = '\t')
colnames(down) = c('Gene', 'Time', 'DE')

#-------------splitting data by DE-----------------------------
up_true_up = subset(up, DE == 'up')
up_true_down = subset(up, DE == 'down')

down_true_up = subset(down, DE == 'up')
down_true_down = subset(down, DE == 'down')

final_up = rbind(up_true_up, down_true_up)
final_up %>% distinct() -> final_up_unique

final_down = rbind(up_true_down, down_true_down)
final_down %>% distinct() -> final_down_unique

#--------------splitting UP DEG by timecourse-----------------------
str(final_up_unique)

final_up_05 = subset(final_up_unique, Time == 0.5)
final_up_3 = subset(final_up_unique, Time == 3.0)
final_up_6 = subset(final_up_unique, Time == 6.0)
final_up_12 = subset(final_up_unique, Time == 12.0)
final_up_24 = subset(final_up_unique, Time == 24.0)
final_up_48 = subset(final_up_unique, Time == 48.0)

#-------------splitting DOWN DEG by timecourse------------------------

str(final_down_unique)

final_down_05 = subset(final_down_unique, Time == 0.5)
final_down_3 = subset(final_down_unique, Time == 3.0)
final_down_6 = subset(final_down_unique, Time == 6.0)
final_down_12 = subset(final_down_unique, Time == 12.0)
final_down_24 = subset(final_down_unique, Time == 24.0)
final_down_48 = subset(final_down_unique, Time == 48.0)

#------intersecting with repeats nearby all genes---------

all_repeats = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/all_repeats_gene_annotation_filtered_dist.tsv',
                   stringsAsFactors = F, header = F, sep = '\t')
str(all_repeats)
names(all_repeats)[names(all_repeats) == 'V20'] <- 'Gene'
all_repeats$Repeats = sapply(all_repeats$V3, function(x) {strsplit(x, '_')[[1]][1]})
all_repeats$Repeats[all_repeats$Repeats == 'Low'] = 'Low_complexity'
all_repeats$Repeats[all_repeats$Repeats == 'Simple'] = 'Simple_repeat'
all_repeats = subset(all_repeats, Repeats != 'Low_complexity' & Repeats != 'Simple_repeat'
                     & Repeats != 'Unknown')

#--------all repeats-------------

freq_rep_all = as.data.frame(table(all_repeats$Repeats))
freq_rep_all$Var1 = as.character(freq_rep_all$Var1)

# intersecting with all repeats from repeatcraft

all_repeats_rtpcrft = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/hvulg_rptcrft_merged.rmerge_renamed_sorted.gff',
                   stringsAsFactors = F, header = F, sep = '\t')
#head(View(all_repeats))
all_repeats_rtpcrft$Repeats = sapply(all_repeats_rtpcrft$V3, function(x) {strsplit(x, '_')[[1]][1]})
all_repeats_rtpcrft$Repeats[all_repeats_rtpcrft$Repeats == 'Low'] = 'Low_complexity'
all_repeats_rtpcrft$Repeats[all_repeats_rtpcrft$Repeats == 'Simple'] = 'Simple_repeat'

#counting frequency
freq_rep_all_rpt = as.data.frame(table(all_repeats_rtpcrft$Repeats))
freq_rep_all_rpt$Var1 = as.character(freq_rep_all_rpt$Var1)
freq_rep_all_rpt = subset(freq_rep_all_rpt, Freq > 200)
freq_rep_all_rpt = subset(freq_rep_all_rpt, Var1 != 'Unknown')


#----------0.5 hours regeneration-----------------
all_repeats_up_05 <- merge(final_up_05, all_repeats, by='Gene')
rep_up_freq_05 = as.data.frame(table(all_repeats_up_05$Repeats))
rep_up_freq_05$Var1 = as.character(rep_up_freq_05$Var1)
rep_up_freq_05_merged = merge(rep_up_freq_05, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_05_merged$Repeats_norm = rep_up_freq_05_merged$Freq.x/rep_up_freq_05_merged$Freq.y/length(final_up_05$Gene)


all_repeats_down_05 <- merge(final_down_05, all_repeats, by='Gene')
rep_down_freq_05 = as.data.frame(table(all_repeats_down_05$Repeats))
rep_down_freq_05$Var1 = as.character(rep_down_freq_05$Var1)
rep_down_05_freq_merged = merge(rep_down_freq_05, freq_rep_all_rpt, by = 'Var1')
rep_down_05_freq_merged$Repeats_norm = rep_down_05_freq_merged$Freq.x/rep_down_05_freq_merged$Freq.y/length(final_down_05$Gene)

#----------3 hours regeneration-----------------
all_repeats_up_3 <- merge(final_up_3, all_repeats, by='Gene')
rep_up_freq_3 = as.data.frame(table(all_repeats_up_3$Repeats))
rep_up_freq_3$Var1 = as.character(rep_up_freq_3)
rep_up_freq_3_merged = merge(rep_up_freq_3, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_3_merged$Repeats_norm = rep_up_freq_3_merged$Freq.x/rep_up_freq_3_merged$Freq.y/length(final_up_3$Gene)

all_repeats_down_3 <- merge(final_down_3, all_repeats, by='Gene')
rep_down_freq_3 = as.data.frame(table(all_repeats_down_3$Repeats))
rep_down_freq_3$Var1 = as.character(rep_down_freq_3$Var1)
rep_down_freq_3_merged = merge(rep_down_freq_3, freq_rep_all_rpt, by = 'Var1')
rep_down_freq_3_merged$Repeats_norm = rep_down_freq_3_merged$Freq.x/rep_down_freq_3_merged$Freq.y/length(final_down_3$Gene)

#----------6 hours of regeneration------------
all_repeats_up_6 <- merge(final_up_6, all_repeats, by='Gene')
rep_up_freq_6 = as.data.frame(table(all_repeats_up_6$Repeats))
rep_up_freq_6$Var1 = as.character(rep_up_freq_6$Var1)
rep_up_freq_6_merged = merge(rep_up_freq_6, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_6_merged$Repeats_norm = rep_up_freq_6_merged$Freq.x/rep_up_freq_6_merged$Freq.y/length(final_up_6$Gene)

all_repeats_down_6 <- merge(final_down_6, all_repeats, by='Gene')
rep_down_freq_6 = as.data.frame(table(all_repeats_down_6$Repeats))
rep_down_freq_6$Var1 = as.character(rep_down_freq_6$Var1)
rep_down_freq_6_merged = merge(rep_down_freq_6, freq_rep_all_rpt, by = 'Var1')
rep_down_freq_6_merged$Repeats_norm = rep_down_freq_6_merged$Freq.x/rep_down_freq_6_merged$Freq.y/length(final_down_6$Gene)

#----------12 hours regeneration-----------------
all_repeats_up_12 <- merge(final_up_12, all_repeats, by='Gene')
rep_up_freq_12 = as.data.frame(table(all_repeats_up_12$Repeats))
rep_up_freq_12$Var1 = as.character(rep_up_freq_12$Var1)
rep_up_freq_12_merged = merge(rep_up_freq_12, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_12_merged$Repeats_norm = rep_up_freq_12_merged$Freq.x/rep_up_freq_12_merged$Freq.y/length(final_up_12$Gene)


all_repeats_down_12 <- merge(final_down_12, all_repeats, by='Gene')
rep_down_freq_12 = as.data.frame(table(all_repeats_down_12$Repeats))
rep_down_freq_12$Var1 = as.character(rep_down_freq_12$Var1)
rep_down_freq_12_merged = merge(rep_down_freq_12, freq_rep_all_rpt, by = 'Var1')
rep_down_freq_12_merged$Repeats_norm = rep_down_freq_12_merged$Freq.x/rep_down_freq_12_merged$Freq.y/length(final_down_12$Gene)

#----------24 hours regeneration-----------------
all_repeats_up_24 <- merge(final_up_24, all_repeats, by='Gene')
rep_up_freq_24 = as.data.frame(table(all_repeats_up_24$Repeats))
rep_up_freq_24$Var1 = as.character(rep_up_freq_24$Var1)
rep_up_freq_24_merged = merge(rep_up_freq_24, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_24_merged$Repeats_norm = rep_up_freq_24_merged$Freq.x/rep_up_freq_24_merged$Freq.y/length(final_up_24$Gene)


all_repeats_down_24 <- merge(final_down_24, all_repeats, by='Gene')
rep_down_freq_24 = as.data.frame(table(all_repeats_down_24$Repeats))
rep_down_freq_24$Var1 = as.character(rep_down_freq_24$Var1)
rep_down_freq_24_merged = merge(rep_down_freq_24, freq_rep_all_rpt, by = 'Var1')
rep_down_freq_24_merged$Repeats_norm = rep_down_freq_24_merged$Freq.x/rep_down_freq_24_merged$Freq.y/length(final_down_24$Gene)

#----------48 hours regeneration-----------------
all_repeats_up_48 <- merge(final_up_48, all_repeats, by='Gene')
rep_up_freq_48 = as.data.frame(table(all_repeats_up_48$Repeats))
rep_up_freq_48$Var1 = as.character(rep_up_freq_48$Var1)
rep_up_freq_48_merged = merge(rep_up_freq_48, freq_rep_all_rpt, by = 'Var1')
rep_up_freq_48_merged$Repeats_norm = rep_up_freq_48_merged$Freq.x/rep_up_freq_48_merged$Freq.y/length(final_up_48$Gene)


all_repeats_down_48 <- merge(final_down_48, all_repeats, by='Gene')
rep_down_freq_48 = as.data.frame(table(all_repeats_down_48$Repeats))
rep_down_freq_48$Var1 = as.character(rep_down_freq_48$Var1)
rep_down_freq_48_merged = merge(rep_down_freq_48, freq_rep_all_rpt, by = 'Var1')
rep_down_freq_48_merged$Repeats_norm = rep_down_freq_48_merged$Freq.x/rep_down_freq_48_merged$Freq.y/length(final_down_48$Gene)

#----------joining dataframes to plot heatmap---------

up_down_05 = full_join(rep_up_freq_05_merged[, c(1,4)], rep_down_05_freq_merged[, c(1,4)],
              by = 'Var1', suffix = c('.05up', '.05down'))

up_down_3 = full_join(rep_up_freq_3_merged[, c(1,4)], rep_down_freq_3_merged[, c(1,4)],
                       by = 'Var1', suffix = c('.3up', '.3down'))

up_down_6 = full_join(rep_up_freq_6_merged[, c(1,4)], rep_down_freq_6_merged[, c(1,4)],
                      by = 'Var1', suffix = c('.6up', '.6down'))

up_down_12 = full_join(rep_up_freq_12_merged[, c(1,4)], rep_down_freq_12_merged[, c(1,4)],
                      by = 'Var1', suffix = c('.12up', '.12down'))

up_down_24 = full_join(rep_up_freq_24_merged[, c(1,4)], rep_down_freq_24_merged[, c(1,4)],
                      by = 'Var1', suffix = c('.24up', '.24down'))

up_down_48 = full_join(rep_up_freq_48_merged[, c(1,4)], rep_down_freq_48_merged[, c(1,4)],
                       by = 'Var1', suffix = c('.48up', '.48down'))

up_down_05_3 = full_join(up_down_05, up_down_3, by = 'Var1',
                         suffix = c('', ''))

up_down_6_12 = full_join(up_down_6, up_down_12, by = 'Var1',
                         suffix = c('', ''))

up_down_24_48 = full_join(up_down_24, up_down_48, by = 'Var1',
                         suffix = c('', ''))

up_down_05_3_6_12 = full_join(up_down_05_3, up_down_6_12, by = 'Var1',
                              suffix = c('', ''))

up_down_05_3_6_12_24_48 = full_join(up_down_05_3_6_12, up_down_24_48, by = 'Var1',
                                    suffix = c('', ''))

up_down_05_3_6_12_24_48[is.na(up_down_05_3_6_12_24_48)] = 0

colnames(up_down_05_3_6_12_24_48) = c('Repeat', 'Up 0.5h', 'Down 0.5h',
                                      'Up 3h', 'Down 3h',
                                      'Up 6h', 'Down 6h',
                                      'Up 12h', 'Down 12h',
                                      'Up 24h', 'Down 24h',
                                      'Up 48h', 'Down 48h')

rownames(up_down_05_3_6_12_24_48) = up_down_05_3_6_12_24_48$Repeat

up_down_05_3_6_12_24_48_melt <- gather(up_down_05_3_6_12_24_48[, -1], Time, Value)
up_down_05_3_6_12_24_48_melt$'Regulation type' = sapply(up_down_05_3_6_12_24_48_melt$Time, function(x) {strsplit(x, ' ')[[1]][1]})
up_down_05_3_6_12_24_48_melt$Time = as.factor(up_down_05_3_6_12_24_48_melt$Time)
str(up_down_05_3_6_12_24_48_melt)
up_down_05_3_6_12_24_48_melt$Time <- ordered(up_down_05_3_6_12_24_48_melt$Time,
                                             levels = c('Up 0.5h', 'Up 3h', 'Up 6h', 'Up 12h', 'Up 24h', 'Up 48h',
                                                        'Down 0.5h', 'Down 3h', 'Down 6h', 'Down 12h', 'Down 24h', 'Down 48h'))



boxplot_time = ggplot(up_down_05_3_6_12_24_48_melt, aes(Time, Value, color = `Regulation type`)) +
  geom_boxplot(fill = 'gray88') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 60, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size =18,color = "black")) +
  labs(title = 'Percentage of repeats per differentially expressed gene in each stage',
       x = 'Regeneration stages', y = 'Percentage') +
  scale_y_continuous(breaks=seq(0, 0.0004, 0.0001), limits=c(0, 0.00045)) +
  scale_color_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2]))

ggsave(plot = boxplot_time, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/boxplot_deg.tiff',
       width = 10, height = 10, device = 'tiff')

up_down_05_3_6_12_24_48 = up_down_05_3_6_12_24_48[, c(1,2,4,6,8,10,12,3,5,7,9,11,13)]
heatmap_up_down = pheatmap(up_down_05_3_6_12_24_48[, -1],
                           cellwidth = 20, cellheight = 10, angle_col = 45,
                           scale = 'row',
                           cluster_cols = F,
                           color = colorRampPalette(c(brewer.pal(9, name = "Set1")[2],
                                                      'white',
                                                      brewer.pal(9, name = "Set1")[1]))(50),
                           fontsize = 14, fontsize_row = 9, fontsize_col = 14,
                           clustering_distance_rows = "manhattan",
                          # clustering_distance_cols = "correlation"
                          legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.15),
                          legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "Repeat number (Z-score)\n")
)

ggsave(plot = heatmap_up_down,
       filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/heatmap_up_down_stages.tiff',
       width = 10, height = 10, device = 'tiff')
ggsave(plot = heatmap_up_down,
       filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/heatmap_up_down_stages.jpeg',
       width = 10, height = 10, device = 'jpeg')

up_down_05_3_6_12_24_48_t = as.data.frame(t(up_down_05_3_6_12_24_48[, -1]))
up_down_05_3_6_12_24_48_t$DE = sapply(rownames(up_down_05_3_6_12_24_48_t),
                                           function(x) {strsplit(x, ' ')[[1]][1]})
str(up_down_05_3_6_12_24_48_t)
pca_res <- prcomp(up_down_05_3_6_12_24_48_t[, -60], scale. = TRUE)

pca = autoplot(pca_res, data = up_down_05_3_6_12_24_48_t,
               label.size = 5, label = T,
               shape = FALSE, colour = 'DE', label.hjust = 0.5,
               label.vjust = 0.5
               
) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.position = 'none'
        ) +
  labs(title = 'PCA plot of regeneration stages',
       x = 'PC1 (29.23%)', y = 'PC2 (17.97%)') + 
  scale_color_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2]))


ggsave(plot = pca, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pca_deg_rep.tiff',
       width = 10, height = 10, device = 'tiff')
ggsave(plot = pca, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pca_deg_rep.jpeg',
       width = 10, height = 10, device = 'jpeg')

# positive and negative contribution
# PC1 positive
rotation_val = as.data.frame(pca_res$rotation)
rotation_pos_pc1 = subset(rotation_val, PC1 > 0)
rotation_pos_pc1 = rotation_pos_pc1[order(-rotation_pos_pc1$PC1),]
pc1_rotation_pos_pc1_top5 = data.frame('PC1' = rotation_pos_pc1$PC1[1:5],
                                        'Repeat' = rownames(rotation_pos_pc1)[1:5])
rotation_neg_pc1 = subset(rotation_val, PC1 < 0)
rotation_neg_pc1 = rotation_neg_pc1[order(rotation_neg_pc1$PC1),]
pc1_rotation_neg_pc1_top5 = data.frame('PC1' = rotation_neg_pc1$PC1[1:5],
                                        'Repeat' = rownames(rotation_neg_pc1)[1:5])
pc1_pos_neg_top5 = rbind(pc1_rotation_pos_pc1_top5, pc1_rotation_neg_pc1_top5)
pc1_pos_neg_top5$Value = rep(c('Positive', 'Negative'), each = 5)

pc1 = ggplot(pc1_pos_neg_top5, aes(reorder(Repeat, -PC1), PC1))+
  geom_bar(stat = 'identity', position = 'dodge', fill = brewer.pal(8, 'Pastel2')[7], color = 'gray55') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 60, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size =18,color = "black")) +
  labs(title = 'Contributors to PC1 divergence',
       x = '', y = 'PC1') +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), limits=c(-0.5, 0.5))
ggsave(plot = pc1, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pc1_div.tiff',
       width = 10, height = 10, device = 'tiff')

pc1_avg <- pc1_pos_neg_top5 %>%
  summarize(avg = mean(PC1, na.rm = T)) %>%
  pull(avg)

str(pc1_pos_neg_top5)

pc1_pretty = ggplot(pc1_pos_neg_top5, aes(reorder(Repeat, -PC1), PC1, color= Value)) +
  geom_segment(aes(x = reorder(Repeat, -PC1), xend = Repeat,
                   y = 0, yend = PC1),
               size = 0.7, alpha = 0.8) +
  geom_hline(aes(yintercept = 0), color = "gray55", size = 0.5) +
  geom_point(aes(reorder(Repeat, -PC1)), size = 5, alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18, color = 'black'),
        axis.text.x = element_text(size = 18,color = 'black'),
        axis.text.y = element_text(size = 18, color = 'black'),
  ) +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), limits=c(-0.4, 0.4)) +
  labs(x = NULL, y = "PC1") +
  coord_flip() +
  scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[1], brewer.pal(11, 'RdYlGn')[11]))

ggsave(plot = pc1_pretty, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pc1_div_pretty.tiff',
       width = 10, height = 10, device = 'tiff')

# PC2 pos and neg
rotation_val = as.data.frame(pca_res$rotation)
rotation_pos_pc2 = subset(rotation_val, PC2 > 0)
rotation_pos_pc2 = rotation_pos_pc2[order(-rotation_pos_pc2$PC2),]
pc2_rotation_pos_top5 = data.frame('PC2' = rotation_pos_pc2$PC2[1:5],
                                       'Repeat' = rownames(rotation_pos_pc2)[1:5])
rotation_neg_pc2 = subset(rotation_val, PC2 < 0)
rotation_neg_pc2 = rotation_neg_pc2[order(rotation_neg_pc2$PC2),]
pc2_rotation_neg_top5 = data.frame('PC2' = rotation_neg_pc2$PC2[1:5],
                                       'Repeat' = rownames(rotation_neg_pc2)[1:5])
pc2_pos_neg_top5 = rbind(pc2_rotation_pos_top5, pc2_rotation_neg_top5)
pc2_pos_neg_top5$Value = rep(c('Positive', 'Negative'), each = 5)

pc2 = ggplot(pc2_pos_neg_top5, aes(reorder(Repeat, -PC2), PC2))+
  geom_bar(stat = 'identity', position = 'dodge', fill = brewer.pal(8, 'Pastel1')[7], color = 'gray55') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 60, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size =18,color = "black")) +
  labs(title = 'Contributors to PC2 divergence',
       x = '', y = 'PC2') +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), limits=c(-0.5, 0.5))

ggsave(plot = pc2, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pc2_div.tiff',
       width = 10, height = 10, device = 'tiff')


pc2_avg <- pc2_pos_neg_top5 %>%
  summarize(avg = mean(PC2, na.rm = T)) %>%
  pull(avg)

pc2_pretty = ggplot(pc2_pos_neg_top5, aes(reorder(Repeat, -PC2), PC2, color= Value)) +
  geom_segment(aes(x = reorder(Repeat, -PC2), xend = Repeat,
                 y = 0, yend = PC2),
             size = 0.7, alpha = 0.8) +
  geom_hline(aes(yintercept = 0), color = "gray55", size = 0.5) +
  geom_point(aes(reorder(Repeat, -PC2)), size = 5, alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18, color = 'black'),
        axis.text.x = element_text(size = 18, color = 'black'),
        axis.text.y = element_text(size = 18, color = 'black'),
        ) +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), limits=c(-0.4, 0.4)) +
  labs(x = NULL, y = "PC2") +
  coord_flip() +
  scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[1], brewer.pal(11, 'RdYlGn')[11]))

ggsave(plot = pc2_pretty, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats_deg/pc2_div_pretty.tiff',
       width = 10, height = 10, device = 'tiff')
