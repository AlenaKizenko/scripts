
library(ggplot2)
library(plyr)
library(RColorBrewer)


#-----------------------------hydra vulgaris repeats------------------------------------------------------------

df_blast_out_hvir = read.table('/home/alena/Documents/PhD/green_hydra_repeat_ann/blast_CR1_proteins/hvir_trans_blast_swiss', sep = '\t')
colnames(df_blast_out_hvir) = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
df_blast_out_hvir_filt = subset(df_blast_out_hvir, evalue < 0.01)
protein_hits_hvir = sapply(df_blast_out_hvir_filt$sseqid, FUN = function(x) {
  strsplit(as.character(x), "\\.")[[1]][1]
})

protein_hits_hvir = unique(unlist(protein_hits_hvir))
#write.csv(protein_hits_hvir, '/home/alena/Documents/PhD/green_hydra_repeat_ann/blast_CR1_proteins/protein_hits.csv')

#performing UniProt search

df_unihits_hvir = read.table('/home/alena/Documents/PhD/green_hydra_repeat_ann/blast_CR1_proteins/uniprot-hits_hvir.tab', sep = '\t', header = T)


#-------------------hydra oligactis repeats---------------------------------------------

df_blast_out_holig = read.table('/home/alena/Documents/PhD/hydra_oligactis/blast_CR1_proteins/holig_trans_blast_swiss', sep = '\t')
colnames(df_blast_out_holig) = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
df_blast_out_holig_filt = subset(df_blast_out_holig, evalue < 0.01)

protein_hits_holig = sapply(df_blast_out_holig_filt$sseqid, FUN = function(x) {
  strsplit(as.character(x), "\\.")[[1]][1]
  })

protein_hits_holig = unique(unlist(protein_hits_holig))

#write.csv(protein_hits_holig, '/home/alena/Documents/PhD/hydra_oligactis/blast_CR1_proteins/protein_hits_holig.csv')

#performing UniProt search

df_unihits_holig = read.csv('/home/alena/Documents/PhD/hydra_oligactis/blast_CR1_proteins/uniprot-hits_holig.tab', sep = '\t', header = T)

#-----------------------hydra vulgaris chr scale assembly repeats----------------------------------------------

df_blast_out_hvulg = read.table('/home/alena/Documents/PhD/hydra_vulgaris_chrom/blast_CR1_proteins/hvulg_trans_blast_swiss', sep = '\t')
colnames(df_blast_out_hvulg) = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
df_blast_out_hvulg_filt = subset(df_blast_out_hvulg, evalue < 0.01)
protein_hits_hvulg = sapply(df_blast_out_hvulg_filt$sseqid, FUN = function(x) {
  strsplit(as.character(x), "\\.")[[1]][1]
})

protein_hits_hvulg = unique(unlist(protein_hits_hvulg))
#write.csv(protein_hits_hvulg, '/home/alena/Documents/PhD/hydra_vulgaris_chrom/blast_CR1_proteins/protein_hits_hvulgchr.csv')

#performing UniProt search

df_unihits_hvulgchr = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/blast_CR1_proteins/uniprot-hits_hvulg.tab', sep = '\t', header = T)

# plotting GO molecular functions
# plotting Hydra viridissima proteins

protein_mf_hvir = unlist(sapply(df_unihits_hvir$Gene.ontology..molecular.function., FUN = function(x) {
  strsplit(as.character(x), "; ")
}))


df_mf_hvir = data.frame('Molecular_function' = c(protein_mf_hvir), 'Species' = rep('Hydra_viridissima', each = length(c(protein_mf_hvir))))
mf_freq_hvir = count(df_mf_hvir$Molecular_function)
mf_freq_hvir = mf_freq_hvir[order(-mf_freq_hvir$freq),]
mf_freq_hvir$x = unlist(sapply(mf_freq_hvir$x, FUN = function(x) {
  strsplit(as.character(x), " \\[")[[1]][1]
}))
nrow(mf_freq_hvir)/5

hvir_GO = ggplot(mf_freq_hvir[1:2,], aes(reorder(x, freq), freq)) +
  geom_histogram(stat = 'identity', fill = brewer.pal(3, 'Dark2')[2], color = 'gray8') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black')
  ) +
  labs(x = 'GO molecular function', y = 'Number of proteins') +
  coord_flip()

# plotting Hydra oligactis proteins

protein_mf_holig = unlist(sapply(df_unihits_holig$Gene.ontology..molecular.function., FUN = function(x) {
  strsplit(as.character(x), "; ")
}))


df_mf_holig = data.frame('Molecular_function' = c(protein_mf_holig), 'Species' = rep('Hydra_oligactis', each = length(c(protein_mf_holig))))
mf_freq_holig = count(df_mf_holig$Molecular_function)
mf_freq_holig = mf_freq_holig[order(-mf_freq_holig$freq),]
mf_freq_holig$x = unlist(sapply(mf_freq_holig$x, FUN = function(x) {
  strsplit(as.character(x), " \\[")[[1]][1]
}))
nrow(mf_freq_holig)/5

holig_GO = ggplot(mf_freq_holig[1:32,], aes(reorder(x, freq), freq)) +
  geom_histogram(stat = 'identity', fill = brewer.pal(3, 'Dark2')[2], color = 'gray8') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black')
  ) +
  labs(x = 'GO molecular function', y = 'Number of proteins') +
  coord_flip()


# plotting Hydra vulgaris proteins

protein_mf_hvulg = unlist(sapply(df_unihits_hvulgchr$Gene.ontology..molecular.function., FUN = function(x) {
  strsplit(as.character(x), "; ")
}))


df_mf_hvulg = data.frame('Molecular_function' = c(protein_mf_hvulg), 'Species' = rep('Hydra_vulgaris', each = length(c(protein_mf_hvulg))))
mf_freq_hvulg = count(df_mf_hvulg$Molecular_function)
mf_freq_hvulg = mf_freq_hvulg[order(-mf_freq_hvulg$freq),]
mf_freq_hvulg$x = unlist(sapply(mf_freq_hvulg$x, FUN = function(x) {
  strsplit(as.character(x), " \\[")[[1]][1]
}))

nrow(mf_freq_hvulg)/5

hvulg_GO = ggplot(mf_freq_hvulg[1:28,], aes(reorder(x, freq), freq)) +
  geom_histogram(stat = 'identity', fill = brewer.pal(3, 'Dark2')[3]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black')
  ) +
  labs(x = 'GO molecular function', y = 'Number of proteins') +
  coord_flip()

#-------------------------Venn diagram----------------------------------------------

library(VennDiagram)
myCol <- brewer.pal(3, "Dark2")

hvir_GO_venn = c(mf_freq_hvir$x)
holig_GO_venn = c(mf_freq_holig$x)
hvulg_GO_venn = c(mf_freq_hvulg$x)


venn.plot = venn.diagram(x = list(holig_GO_venn, hvulg_GO_venn, hvir_GO_venn),
             category.names = c('Hydra oligactis', 'Hydra vulgaris', 'Hydra viridissima'),
             filename = NULL,
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             output = T,
             # Numbers
             cex = 1,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
            # cat.pos = c(-25, 27, 160),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)
dev.off()
grid.draw(venn.plot)

myCol <- c("#D95F02", "#1B9E77", "#7570B3")

entries_hvir = as.character(df_unihits_hvir$Entry)
entries_holig = as.character(df_unihits_holig$Entry)
entries_hvulg = as.character(df_unihits_hvulgchr$Entry)

venn.plot = venn.diagram(x = list(entries_hvir, entries_holig, entries_hvulg),
                         category.names = c('Hydra viridissima', 'Hydra oligactis', 'Hydra vulgaris'),
                         filename = NULL,
                         lwd = 2,
                         lty = 'blank',
                         fill = myCol,
                         output = T,
                         # Numbers
                         cex = 1.5,
                         fontface = "bold",
                         fontfamily = "sans",
                         
                         # Set names
                         cat.cex = 1.5,
                         cat.fontface = "bold",
                         cat.default.pos = "outer",
                         # cat.pos = c(-25, 27, 160),
                         cat.dist = c(0.055, 0.055, 0.085),
                         cat.fontfamily = "sans",
                         rotation = 1)
dev.off()
grid.draw(venn.plot)


setdiff(intersect(entries_holig, entries_hvulg), intersect(entries_hvir, entries_holig))

intersect(entries_hvir, entries_holig)
intersect(entries_hvir, entries_hvulg)


#----------------------------intersection of H. oligactis, H. viridissima, H. vulgarischr proteins--------------------------

cr1_intersection = read.csv('/home/alena/Documents/PhD/CR1_proteins/H_oli_vir_vulchr_intersection_proteins.tab', sep = '\t', header = T)


protein_mf_all = unlist(sapply(cr1_intersection$Gene.ontology..molecular.function., FUN = function(x) {
  strsplit(as.character(x), "; ")
}))


df_mf_all = data.frame('Molecular_function' = c(protein_mf_all))
mf_freq_all = count(df_mf_all$Molecular_function)
mf_freq_all = mf_freq_all[order(-mf_freq_all$freq),]
mf_freq_all$x = unlist(sapply(mf_freq_all$x, FUN = function(x) {
  strsplit(as.character(x), " \\[")[[1]][1]
}))


intersect_GO = ggplot(mf_freq_all, aes(reorder(x, freq), freq)) +
  geom_histogram(stat = 'identity', fill = brewer.pal(8, 'Accent')[5]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16)
  ) +
  labs(x = 'GO molecular function', y = 'Number of proteins') +
  coord_flip()



#----------------------------intersection of H. oligactis and H. vulgarischr proteins--------------------------

cr1_int_olivul = read.csv('/home/alena/Documents/PhD/CR1_proteins/H_oli_vulchr_intersection_proteins.tab', sep = '\t', header = T)


protein_mf_all = unlist(sapply(cr1_int_olivul$Gene.ontology..molecular.function., FUN = function(x) {
  strsplit(as.character(x), "; ")
}))


df_mf_all = data.frame('Molecular_function' = c(protein_mf_all))
mf_freq_all = count(df_mf_all$Molecular_function)
mf_freq_all = mf_freq_all[order(-mf_freq_all$freq),]
mf_freq_all$x = unlist(sapply(mf_freq_all$x, FUN = function(x) {
  strsplit(as.character(x), " \\[")[[1]][1]
}))


intersect_GO = ggplot(mf_freq_all, aes(reorder(x, freq), freq)) +
  geom_histogram(stat = 'identity', fill = brewer.pal(8, 'Accent')[7]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16)
  ) +
  labs(x = 'GO molecular function', y = 'Number of proteins') +
  coord_flip()
