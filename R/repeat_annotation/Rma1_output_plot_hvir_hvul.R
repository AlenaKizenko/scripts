library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(viridis)


#---------------Hydra viridissima-----------------------------
lmig_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_viridissima/round_1/repeatmasker_parse_result/hvir_divsum.txt',
                   header = T, sep = '\t', skip = 4)

lmig_hvir_families = lmig_hvir[c(2:2340),]
cr1_hvir = subset(lmig_hvir_families, Class == 'LINE/CR1')
cr1_hvir$Kimura. = as.numeric(cr1_hvir$Kimura.)

start_num = which(lmig_hvir$Class == 'Coverage for each repeat class and divergence (Kimura)')
lmig_hvir = as.data.frame(lmig_hvir[-c(1:start_num),1])
colnames(lmig_hvir) = 'V1'
names = as.vector(strsplit(lmig_hvir$V1[1], split = ' ')[[1]])
lmig_hvir = separate(lmig_hvir, V1,
                into = as.vector(names),
                sep = ' ')

lmig_hvir = as.data.frame(sapply(lmig_hvir[-1,], as.numeric))
lmig_hvir[, -1] = as.data.frame(sapply(lmig_hvir[, -1], FUN = function(x) {x/284000000*100}))
lm_hvir <- melt(lmig_hvir, id.vars = 0:1)
lm_hvir_dna = lm_hvir[grepl('DNA', lm_hvir$variable),]
lm_hvir$variable = sapply(as.character(lm_hvir$variable), function(x) {
  strsplit(x, split = '/', fixed = T)[[1]][1]
  })


#---------------Hydra vulgaris-----------------------------
lmig_hvul = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/rma1_replanscape/hvul_divsum.txt',
                   header = T, sep = '\t', skip = 4)
start_num = which(lmig_hvul$Class == 'Coverage for each repeat class and divergence (Kimura)')
lmig_hvul = as.data.frame(lmig_hvul[-c(1:start_num),1])
colnames(lmig_hvul) = 'V1'
names = as.vector(strsplit(lmig_hvul$V1[1], split = ' ')[[1]])
lmig_hvul = separate(lmig_hvul, V1,
                into = as.vector(names),
                sep = ' ')

lmig_hvul = as.data.frame(sapply(lmig_hvul[-1,], as.numeric))
lmig_hvul[, -1] = as.data.frame(sapply(lmig_hvul[, -1], FUN = function(x) {x/1250000000*100}))
lm_hvul<- melt(lmig_hvul, id.vars = 0:1)
lm_hvul_dna = lm_hvul[grepl('DNA', lm_hvul$variable),]
lm_hvul$variable = sapply(as.character(lm_hvul$variable), function(x) {
  strsplit(x, split = '/', fixed = T)[[1]][1]
  })

#--------merging dataframes----------

lm_hvir = lm_hvir[, -1]
lm_hvir = rbind(lm_hvir, c('No repeats', 100-sum(lm_hvir$value)))
lm_hvir = aggregate(as.numeric(lm_hvir$value), by=list(Category=lm_hvir$variable), FUN=sum)
lm_hvir$Species = rep('Hydra viridissima', length(lm_hvir$Category))

lm_hvul = lm_hvul[, -1]
lm_hvul = rbind(lm_hvul, c('No repeats', 100-sum(lm_hvul$value)))
lm_hvul = aggregate(as.numeric(lm_hvul$value), by=list(Category=lm_hvul$variable), FUN=sum)
lm_hvul$Species = rep('Hydra vulgaris', length(lm_hvul$Category))

lm = rbind(lm_hvir, lm_hvul)
str(lm)
colnames(lm) = c('variable', 'value', 'Species')
lm_subset = subset(lm, lm$variable != 'Other' & lm$variable != 'snRNA' &
                     lm$variable != 'Satellite' & lm$variable != 'Simple_repeat')
unique(lm_subset$variable)
lm_subset$variable = factor(lm_subset$variable, levels = c('LINE','LTR',
                                              'DNA', 'RC',
                                                'SINE', 'Unknown', 'No repeats'))
lm_subset$Species = as.factor(lm_subset$Species)
lm_subset$value = as.numeric(lm_subset$value)

str(lm)

plot_Rma = ggplot(data=lm_subset, aes(Species, value, fill=variable, group = variable))+
  geom_bar(stat="identity", position="fill", width = 0.5, color = 'black', alpha = 0.8) +
  scale_fill_viridis(discrete = T, name = 'Repeat type') +
  theme_bw() +
  theme(
        axis.text.y = element_text(size=18,color = "black"),
        axis.text.x = element_text(size=18,color = "black", face = 'italic'),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
    
  ) +
  scale_y_continuous(labels=scales::percent)

ggsave(plot = plot_Rma, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/repeat_content_comp.tiff',
       device = 'tiff', width = 10, height = 10)

#---------------DNA repeats---------------------------
lm_hvir_dna = subset(lm_hvir_dna, variable != 'Other/DNA_virus')
lm_hvir_dna = aggregate(as.numeric(lm_hvir_dna$value),
                        by=list(Category=lm_hvir_dna$variable), FUN=sum)
colnames(lm_hvir_dna) = c('variable', 'value')
lm_hvir_dna$Species = rep('Hydra viridissima', length(lm_hvir_dna$variable))

lm_hvul_dna = subset(lm_hvul_dna, variable != 'Other/DNA_virus')
lm_hvul_dna = aggregate(as.numeric(lm_hvul_dna$value),
                        by=list(Category=lm_hvul_dna$variable), FUN=sum)
colnames(lm_hvul_dna) = c('variable', 'value')
lm_hvul_dna$Species = rep('Hydra vulgaris', length(lm_hvul_dna$variable))
lm_dna = rbind(lm_hvir_dna, lm_hvul_dna)

plot_dna = ggplot(lm_dna, aes(reorder(variable, -value), value, fill = Species)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, color = 'black', size = 8),
    axis.text.y = element_text(size=18,color = "black"),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "transparent", color = "transparent")
  ) +
  scale_fill_viridis(alpha = 0.8, discrete = T, option = 'E')

ggsave(plot = plot_dna, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/repeat_content_comp_dna.tiff',
       device = 'tiff', width = 10, height = 10)


