library(dplyr)
library(ggplot2)
library(topGO)
library(RColorBrewer)

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

#------intersecting with genes nearby DNA/Kolobok-T2 repeats---------

kolobok = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_ann_DNA_Merlin_up/DNA_Merlin_gene.tsv',
               stringsAsFactors = F, header = F, sep = '\t')

names(kolobok)[names(kolobok) == 'V19'] <- 'Gene'

#----------0.5 hours regeneration-----------------
kolobok_up_05 <- merge(final_up_05, kolobok, by='Gene')
kolobok_down_05 <- merge(final_down_05, kolobok, by='Gene')

#----------3 hours regeneration-----------------
kolobok_up_3 <- merge(final_up_3, kolobok, by='Gene')
kolobok_down_3 <- merge(final_down_3, kolobok, by='Gene')

#----------6 hours of regeneration------------
kolobok_up_6 <- merge(final_up_6, kolobok, by='Gene')
kolobok_down_6 <- merge(final_down_6, kolobok, by='Gene')

#----------12 hours regeneration-----------------
kolobok_up_12 <- merge(final_up_12, kolobok, by='Gene')
kolobok_down_12 <- merge(final_down_12, kolobok, by='Gene')

#----------24 hours regeneration-----------------
kolobok_up_24 <- merge(final_up_24, kolobok, by='Gene')
kolobok_down_24 <- merge(final_down_24, kolobok, by='Gene')

#----------48 hours regeneration-----------------
kolobok_up_48 <- merge(final_up_48, kolobok, by='Gene')
kolobok_down_48 <- merge(final_down_48, kolobok, by='Gene')

# plotting distribution of repeats

de_norm_rep = data.frame('DE_type' = c('Up 0.5h', 'Down 0.5h',
                                       'Up 3h', 'Down 3h',
                                       'Up 6h', 'Down 6h', 
                                       'Up 12h', 'Down 12h',
                                       'Up 24h', 'Down 24h',
                                       'Up 48h', 'Down 48h'),
                         'Repeats_num_near_deg' = c(
                                                    length(kolobok_up_05$Gene)/length(final_up_05$Gene),
                                                    length(kolobok_down_05$Gene)/length(final_down_05$Gene),
                                                    length(kolobok_up_3$Gene)/length(final_up_3$Gene),
                                                    length(kolobok_down_3$Gene)/length(final_down_3$Gene),
                                                    length(kolobok_up_6$Gene)/length(final_up_6$Gene),
                                                    length(kolobok_down_6$Gene)/length(final_down_6$Gene),
                                                    length(kolobok_up_12$Gene)/length(final_up_12$Gene),
                                                    length(kolobok_down_12$Gene)/length(final_down_12$Gene),
                                                    length(kolobok_up_24$Gene)/length(final_up_24$Gene),
                                                    length(kolobok_down_24$Gene)/length(final_down_24$Gene),
                                                    length(kolobok_up_48$Gene)/length(final_up_48$Gene),
                                                    length(kolobok_down_48$Gene)/length(final_down_48$Gene)
                                                    
                         ))
de_norm_rep$DE_type = as.character(de_norm_rep$DE_type)
de_norm_rep$Reg_type = sapply(de_norm_rep$DE_type, function(x) {strsplit(x, ' ')[[1]][1]})
de_norm_rep$Reg_type[de_norm_rep$Reg_type == 'Down'] = 'Downregulated genes'
de_norm_rep$Reg_type[de_norm_rep$Reg_type == 'Up'] = 'Upregulated genes'

deg = ggplot(de_norm_rep, aes(reorder(DE_type, as.numeric(rownames(de_norm_rep))), Repeats_num_near_deg*100)) +
  geom_rect(data = de_norm_rep,
            aes(fill = Reg_type),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.2) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', fill = brewer.pal(8, 'Pastel1')[7]) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 40, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        strip.text.y = element_text(size = 18,colour = "black"),
        legend.position = 'none'
        # legend.text = element_text(size=18,color = "black"),
        #  legend.title = element_text(size = 18,color = "black")
  ) +
  labs(title = 'Amount of repeats nearby differentially expressed regeneration genes',
       x = '', y = 'Percentage of repeats per one gene',
       fill = 'Repeat type') +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2])) +
  facet_grid(rows = vars(Reg_type))


ggsave(plot = deg, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_ann_DNA_Merlin_up/deg_DNA_Merlin_bars.tiff',
       width = 10, height = 10, device = 'tiff')
