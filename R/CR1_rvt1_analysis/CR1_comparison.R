library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(RColorBrewer)


#-------hydra viridissima 1st round repeats---------------
parse_repeat_r1_hvir = read.csv('/Users/alenakizenko/Documents/PhD/project/hydra_viridissima/round_1/parse_rep_land_rma1/hvir_genome_hm2_250116_renamed.fa.align.landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F)

colnames(parse_repeat_r1_hvir)[4:53] = c(0:49)
parse_repeat_r1_hvir$`Repeat family` = paste(parse_repeat_r1_hvir$Rclass,
                                             parse_repeat_r1_hvir$Rfam, sep="/")
df_repeats_r1_fam_hvir = melt(parse_repeat_r1_hvir[, -c(1,2,3)], id.vars = 'Repeat family')
df_repeats_r1_fam_hvir$perc = c((df_repeats_r1_fam_hvir$value*100)/254000000) # genome size is 254 Mb
df_repeats_r1_fam_hvir = subset(df_repeats_r1_fam_hvir, df_repeats_r1_fam_hvir$value > 0)
df_repeats_r1_fam_cr1_hvir = subset(df_repeats_r1_fam_hvir, `Repeat family` == "LINE/CR1")
df_repeats_r1_fam_cr1_hvir$Species = c(rep('Hydra viridissima', each = length(df_repeats_r1_fam_cr1_hvir$`Repeat family`)))

#-----hydra oligactis 1st round repeats--------------------

parse_repeat_r1_holig = read.csv('/Users/alenakizenko/Documents/PhD/project/hydra_oligactis/rma1/hydra_oligactis.fasta.align.landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F, skip = 1, header = T)

parse_repeat_r1_holig = subset(parse_repeat_r1_holig, parse_repeat_r1_holig$Rclass != 'SINE?')
colnames(parse_repeat_r1_holig)[4:53] = c(0:49)
parse_repeat_r1_holig$`Repeat family` = paste(parse_repeat_r1_holig$Rclass, parse_repeat_r1_holig$Rfam, sep="/")
parse_repeat_r1_holig[,-c(1,2,3,54)] = sapply(parse_repeat_r1_holig[,-c(1,2,3,54)], as.numeric )
df_repeats_r1_fam_holig = melt(parse_repeat_r1_holig[, -c(1,2,3)], id.vars = 'Repeat family')
df_repeats_r1_fam_holig$perc = c((df_repeats_r1_fam_holig$value*100)/1450000000) # genome size is 1450 Mb
df_repeats_r1_fam_holig = subset(df_repeats_r1_fam_holig, df_repeats_r1_fam_holig$value > 0)
df_repeats_r1_fam_cr1_holig = subset(df_repeats_r1_fam_holig, `Repeat family` == "LINE/CR1")
df_repeats_r1_fam_cr1_holig$Species = c(rep('Hydra oligactis', each = length(df_repeats_r1_fam_cr1_holig$`Repeat family`)))

#-----hydra vulgaris 1st round repeats--------------------

parse_repeat_r1_hvul = read.csv('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/rma1_replanscape/Hydra-v2b.fasta.align.landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F, skip = 1, header = T)

colnames(parse_repeat_r1_hvul)[4:53] = c(0:49)
parse_repeat_r1_hvul$`Repeat family` = paste(parse_repeat_r1_hvul$Rclass, parse_repeat_r1_hvul$Rfam, sep="/")
df_repeats_r1_fam_hvul = melt(parse_repeat_r1_hvul[, -c(1,2,3)], id.vars = 'Repeat family')
df_repeats_r1_fam_hvul$perc = c((df_repeats_r1_fam_hvul$value*100)/1250000000) # genome size is 1250 Mb
df_repeats_r1_fam_hvul = subset(df_repeats_r1_fam_hvul, df_repeats_r1_fam_hvul$value > 0)
df_repeats_r1_fam_cr1_hvul = subset(df_repeats_r1_fam_hvul, `Repeat family` == "LINE/CR1")
df_repeats_r1_fam_cr1_hvul$Species = c(rep('Hydra vulgaris', each = length(df_repeats_r1_fam_cr1_hvul$`Repeat family`)))

#--------merging all dataframes------------

df_all_hydras = rbind(df_repeats_r1_fam_cr1_hvir, df_repeats_r1_fam_cr1_holig, df_repeats_r1_fam_cr1_hvul)
str(df_all_hydras)
df_all_hydras$Species = as.factor(df_all_hydras$Species)
df_all_hydras$Species <- ordered(df_all_hydras$Species, levels = c('Hydra vulgaris', 'Hydra oligactis', 'Hydra viridissima'))

sum(df_repeats_r1_fam_cr1_hvul$perc)

plot_all_hydras_CR1 = ggplot(data=df_all_hydras, aes(variable, perc, fill = Species))+
  geom_bar(stat="identity", position = 'stack', alpha = .8, width = 1) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = c(0.67, 0.72),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18, color = 'black'),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1, color = 'black'),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title=element_text(size=18, color = 'black')) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  scale_fill_manual(values = c(brewer.pal(3, 'Dark2')[3], brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1]))

ggsave(plot = plot_all_hydras_CR1,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/plot_all_hydras_CR1.tiff',
       width = 10, height = 10, device = 'tiff')



