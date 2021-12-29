library(seqinr)
library(ggplot2)

lines = read.fasta('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/rvt1_all_fam_renamed.fasta')


lens = lapply(lines, function(x) summary(x)$length)
lens = as.data.frame(unlist(lens))
lens$Repeat = sapply(rownames(lens), function(x) strsplit(x, '_')[[1]][1])
colnames(lens) = c('Length', 'Repeat')

ggplot(lens, aes(Repeat, Length)) +
  geom_violin()
