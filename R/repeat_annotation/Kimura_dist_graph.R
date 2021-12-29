library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(cowplot)


#---------------Hydra vulgaris-----------------------------
lmig_hvul = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/rma1_replanscape/hvul_divsum_copies.txt',
                        header = T, sep = '\t')

cr1 = subset(lmig_hvul, Class == 'LINE/CR1')
cr1$Repeats = paste(cr1$Class, cr1$Repeat, sep = '_')
cr1 = cr1[, -c(1,2,3,4)]
#-----------------

new = read.delim2('/Users/alenakizenko/Documents/PhD/project/CR1_graph_fr/new_repeats.txt', header = F)
new_cr1_hvul = as.data.frame(new[grep('hvul*.', new$V1),])
colnames(new_cr1_hvul) = 'V1'
new_cr1_hvul$Class = sapply(new_cr1_hvul$V1, function(x) {strsplit(x = x, split = '_')[[1]][2]})
new_cr1_hvul$Rnd = sapply(new_cr1_hvul$V1,
                          function(x) {strsplit(x = x, split = '_')[[1]][3]})
new_cr1_hvul$Family = sapply(new_cr1_hvul$V1,
                             function(x) {strsplit(x = x, split = '_')[[1]][4]})

new_cr1_hvul$Repeats = paste(new_cr1_hvul$Class, new_cr1_hvul$Rnd, new_cr1_hvul$Family, sep = '_')
new_cr1_hvul = new_cr1_hvul[, -c(1,2,3,4)]
new_cr1_hvul = as.data.frame(unique(new_cr1_hvul))
colnames(new_cr1_hvul) = 'Repeats'

cr1_5perc_new = merge(cr1, new_cr1_hvul, by = 'Repeats')
cr1_5perc_new$Kimura. = as.numeric(cr1_5perc_new$Kimura.)

new_plot = ggplot(cr1_5perc_new, aes(Kimura.)) +
  geom_histogram(fill = brewer.pal(10, 'Paired')[10], color = 'gray9', binwidth = 1) +
  theme_classic() +
  theme(plot.title = element_text(size = 28),
        axis.text.x = element_text(size=28,color = "black"),
        axis.text.y = element_text(size=28,color = "black"),
        axis.title.y = element_text(size=28,color = "black"),
        axis.title.x = element_text(size=28,color = "black")
        ) +
  ggtitle('New repeats') +
  ylim(0,15) +
  xlim(0, 35) +
  xlab('Kimura substitution level')

#------------------------------------------------------------
old = read.delim2('/Users/alenakizenko/Documents/PhD/project/CR1_graph_fr/old_repeats.txt', header = F)
old_cr1_hvul = as.data.frame(old[grep('hvul*.', old$V1),])
colnames(old_cr1_hvul) = 'V1'
old_cr1_hvul$Class = sapply(old_cr1_hvul$V1, function(x) {strsplit(x = x, split = '_')[[1]][2]})
old_cr1_hvul$Rnd = sapply(old_cr1_hvul$V1,
                          function(x) {strsplit(x = x, split = '_')[[1]][3]})
old_cr1_hvul$Family = sapply(old_cr1_hvul$V1,
                             function(x) {strsplit(x = x, split = '_')[[1]][4]})

old_cr1_hvul$Repeats = paste(old_cr1_hvul$Class, old_cr1_hvul$Rnd, old_cr1_hvul$Family, sep = '_')
old_cr1_hvul = old_cr1_hvul[, -c(1,2,3,4)]
old_cr1_hvul = as.data.frame(unique(old_cr1_hvul))
colnames(old_cr1_hvul) = 'Repeats'

cr1_5perc_old = merge(cr1, old_cr1_hvul, by = 'Repeats')
cr1_5perc_old$Kimura. = as.numeric(cr1_5perc_old$Kimura.)

old_plot = ggplot(cr1_5perc_old, aes(Kimura.)) +
  geom_histogram(fill = brewer.pal(10, 'Paired')[9], color = 'gray9', binwidth = 1) +
  theme_classic() +
  theme(plot.title = element_text(size = 28),
        axis.text.x = element_text(size=28,color = "black"),
        axis.text.y = element_text(size=28,color = "black"),
        axis.title.y = element_text(size=28,color = "black"),
        axis.title.x = element_text(size=28,color = "black")
  ) +
  ggtitle('Old repeats') +
  ylim(0,15) +
  xlim(0, 35) +
  xlab('Kimura substitution level')

new_old_5perc_kimura = plot_grid(new_plot, old_plot)

ggsave(plot = new_old_5perc_kimura,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/plot_5perc_kimura.tiff',
       width = 15, height = 10, device = 'tiff')



rvt1_df = data.frame('Repeat_type' = c('Old repeats', 'New repeats'),
                     'RVT1_number' = c(12, 2456),
                     'TE_number'  = c(9567, 99303))
rvt1_df$RVT1_perc_norm = (rvt1_df$RVT1_number/rvt1_df$TE_number) * 100
str(rvt1_df)

plot_rvt1 = ggplot(rvt1_df, aes(reorder(Repeat_type, as.numeric(rownames(rvt1_df))),
                                RVT1_perc_norm, fill = Repeat_type)) +
  geom_bar(stat = 'identity', width = 0.5, color = 'gray9') +
  theme_bw() +
  theme(
        axis.text.x = element_text(size=28,color = "black"),
        axis.text.y = element_text(size=28,color = "black"),
        axis.title.y = element_text(size=28,color = "black"),
        axis.title.x = element_blank(),
        legend.position = 'none'
  ) +
  ylab('Percentage of reverse transcriptase sequences') +
  scale_fill_manual(values = c(brewer.pal(10, 'Paired')[9], brewer.pal(10, 'Paired')[10]),
                    name = 'Repeat type')

ggsave(plot = plot_rvt1,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/plot_rvt1.tiff',
       width = 15, height = 10, device = 'tiff')


#-----------------alatinidae-----------------

lmig_ala = read.delim2('/Users/alenakizenko/Documents/PhD/project/alatinidae/divsum_alatinidae.txt',
                        header = T, sep = '\t', skip = 4)
start_num = which(lmig_ala$Class == 'Coverage for each repeat class and divergence (Kimura)')
lmig_ala = as.data.frame(lmig_ala[-c(1:start_num),1])
colnames(lmig_ala) = 'V1'
names = as.vector(strsplit(lmig_ala$V1[1], split = ' ')[[1]])
lmig_ala = separate(lmig_ala, V1,
                     into = as.vector(names),
                     sep = ' ')

lmig_ala = as.data.frame(sapply(lmig_ala[-1,], as.numeric))
lmig_ala[, -1] = as.data.frame(sapply(lmig_ala[, -1], FUN = function(x) {x/1250000000*100}))
lmig_ala<- melt(lmig_ala, id.vars = 0:1)
lmig_ala$variable = sapply(as.character(lmig_ala$variable), function(x) {
  strsplit(x, split = '/', fixed = T)[[1]][1]
})


plot_Rma = ggplot(lmig_ala, aes(Div, value, fill = variable)) +
  geom_bar(stat = 'identity')

ggsave(plot = plot_Rma, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/repeat_content_comp.tiff',
       device = 'tiff', width = 10, height = 10)


