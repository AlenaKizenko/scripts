library(ggplot2)


chr_repeat = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/chromosomes_repeat_content.tsv', sep = '\t', stringsAsFactors = F, header = F)

str(chr_repeat)
colnames(chr_repeat) = c('Chromosome', 'Chromosome_length', 'Number_of_repeats')
chr_repeat$Repeats_norm = chr_repeat$Number_of_repeats/chr_repeat$Chromosome_length


chr_rep = ggplot(chr_repeat, aes(x = chr_repeat$Chromosome_length, y = chr_repeat$Repeats_norm*1000)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 20, face = 'bold'),
        axis.text.x = element_text(size=20,color = "black", hjust = 1),
        axis.text.y = element_text(size=20,color = "black"),
        axis.title.x = element_text(size=20,color = "black"),
        axis.title.y = element_text(size=20,color = "black")
        ) +
  labs(title = 'Repeats distribution depending on chromosome length',
       x = 'Chromosome length', y = 'Repeats per 1kb')


ggsave(plot = chr_rep, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/chr_rep_norm.tiff',
       width = 10, height = 10, device = 'tiff')
