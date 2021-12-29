library(ggplot2)
library(RColorBrewer)


#---------
# Green Hydra
green_hydra = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_repeats_near_genes_filtered.tsv',
                       sep = '\t', stringsAsFactors = F, header = F)

colnames(green_hydra)[20] = 'Gene'
colnames(green_hydra)[19] = 'Distance'
cr1 = green_hydra[grepl('LINE/CR1', green_hydra$V3), ]
cr1$Distance = as.numeric(as.character(cr1$Distance))

# Brown Hydra

brown_hydra = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/all_repeats_gene_annotation_filtered_dist.tsv',
                       sep = '\t', stringsAsFactors = F, header = F)

colnames(brown_hydra)[20] = 'Gene'
colnames(brown_hydra)[19] = 'Distance'
cr1_brown = brown_hydra[grepl('LINE/CR1', brown_hydra$V3), ]
cr1_brown$Distance = as.numeric(as.character(cr1_brown$Distance))

# ------ creating dataframe for both hydras---------


cr1_g_b = data.frame('Distance' = c(cr1$Distance, cr1_brown$Distance), 'Hydra species' = c(rep('Hydra viridissima', length(cr1$Distance)),
                                                                                           rep('Hydra vulgaris', length(cr1_brown$Gene))))

str(cr1_g_b)
cr1_g_b$Distance = as.numeric(cr1_g_b$Distance)
cr1_g_b$Hydra.species = as.factor(cr1_g_b$Hydra.species)

cr1_g_b_zero = subset(cr1_g_b, Distance == 0)

dist_plot = ggplot(cr1_g_b, aes(Hydra.species, Distance)) +
  geom_boxplot(fill = c(brewer.pal(11, 'RdBu')[8]),
               color = 'gray8', width = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Distribution of gene-repeat distances',
       x = '', y = 'Distance'
  )

ggsave(plot = dist_plot, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/boxplot_g_b.tiff',
       device = 'tiff', width = 10, height = 10)

counts_zero = as.data.frame(table(cr1_g_b_zero))
counts = as.data.frame(table(cr1_g_b$Hydra.species))
counts_final = cbind(counts_zero, counts)
colnames(counts_final) = c("Distance_zero", "Hydra.species", "Freq_zero", "Distance", "Freq" )
counts_final$Norm = counts_final$Freq_zero/counts_final$Freq

zero_dist = ggplot(counts_final, aes(Distance, Norm)) +
  geom_bar(stat = 'identity', width = 0.5,
           fill = brewer.pal(11, 'RdBu')[8],
           color = 'gray8') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Percentage of repeats inside genes',
       x = '', y = 'Percetage'
  )

ggsave(plot = zero_dist, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/barblot-zero_g_b.tiff',
       device = 'tiff', width = 10, height = 10)

#-------comparing scaffolds distribution-------
scaf_green = as.data.frame(table(cr1$V1))

green_scaffold = ggplot(scaf_green, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = 'identity', width = 1, fill = brewer.pal(11, 'RdBu')[5],
           color = 'gray8') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Distribution of repeats in Hydra viridissima genome',
       x = 'Scaffolds', y = 'Number of repeats'
  )

ggsave(plot = green_scaffold, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/scaff_g.tiff',
       device = 'tiff', width = 10, height = 10)
#brown

scaf_brown = as.data.frame(table(cr1_brown$V1))

brown_scaffold = ggplot(scaf_brown, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = 'identity', width = 1, fill = brewer.pal(11, 'RdBu')[7],
           color = 'gray8') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Distribution of repeats in Hydra vulgaris genome',
       x = 'Scaffolds', y = 'Number of repeats'
  )
ggsave(plot = brown_scaffold, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/scaff_b.tiff',
       device = 'tiff', width = 10, height = 10)


#---------------Hydra viridissima-----------------------------
lmig_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/green_hydra_repeat_ann/round_1/repeatmasker_parse_result/hvir_divsum.txt',
                        header = T, sep = '\t', skip = 4)
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
colnames(lm) = c('variable', 'value', 'Species')
lm_subset = subset(lm, lm$variable != 'Other' & lm$variable != 'snRNA' &
                     lm$variable != 'Satellite' & lm$variable != 'Simple_repeat')
lm_subset$variable = factor(lm_subset$variable, levels = c('LINE','LTR',
                                                           'DNA', 'RC',
                                                           'SINE', 'Unknown', 'No repeats'))
lm_subset$Species = as.factor(lm_subset$Species)
lm_subset$value = as.numeric(lm_subset$value)

plot_Rma = ggplot(data=lm_subset, aes(Species, value, fill=variable, group = variable))+
  geom_bar(stat="identity", position="fill", width = 0.5, color = 'black', alpha = 0.8) +
  scale_fill_viridis(discrete = T, name = 'Repeat type') +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=16,color = "black"),
    axis.text.x = element_text(size=16,color = "black", face = 'italic'),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
    
  ) +
  scale_y_continuous(labels=scales::percent)
plot_Rma