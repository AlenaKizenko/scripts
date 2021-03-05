library(circlize)
library(stats)
library(dendextend)
library(RColorBrewer)
library(igraph)


#-----------------------CR1 1 percent hierarchial clustering--------------------------------------------------------------------------------------------
df_hydras_cr1 = read.table('/home/alena/Documents/PhD/CR1_graphs_all_hydras/h_oli_vir_vulchr_blast_1perc_dist_filtered.tsv',
                           header = T, sep = '\t', stringsAsFactors = F)
df_hydras_cr1 = subset(df_hydras_cr1, Source != Target)
df_hydras_cr1[df_hydras_cr1 == 0] <- 2
g <- graph.data.frame(df_hydras_cr1[, 1:2], directed=F)
E(g)$weight <- df_hydras_cr1$Weight
mat = get.adjacency(g, attr = 'weight', sparse = F)
mat[mat == 0] = 1
mat[mat == 2] = 0

#mat = apply(mat, 2, function(x) {1 - x})
#View(mat)
#diag(mat) = 0
a = hclust(as.dist(mat), method = "average", members = NULL)
dend = as.dendrogram(a)
col_aa_red <- ifelse(grepl("Hviridissima", labels(dend)),
                     brewer.pal(3, 'Dark2')[2],
                     ifelse(grepl("Holigactis", labels(dend)),
                            brewer.pal(3, 'Dark2')[1],
                            brewer.pal(3, 'Dark2')[3]))

dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col")

dend2 %>% set("labels_cex", 20) %>% plot(ylab = 'Distance')
#legend(x = 0.7, y = 0.7, legend = c('Hydra oligactis', 'Hydra viridissima', 'Hydra vulgaris'),
 #      bty = "n", fill=brewer.pal(3, 'Dark2'), border=NA, cex = 1)
title(main = 'Hierarchical clustering of 1% CR1 repeats', cex.main = 2)
dev.off()


#----------------------Penelope 1 percent hierarchical clustering-----------------------------------------------------------------------------------------

df_hydras_penelope1 = read.table('/home/alena/Documents/PhD/Penelope_1perc/h_oli_vul_vir_blast_Penelope1_dist_filtered.tsv',
                                 header = T, sep = '\t')
df_hydras_penelope1 = subset(df_hydras_penelope1, Source != Target)
df_hydras_penelope1[df_hydras_penelope1 == 0] <- 2
g <- graph.data.frame(df_hydras_penelope1[, 1:2], directed=F)
E(g)$weight <- df_hydras_penelope1$Weight
mat = get.adjacency(g, attr = 'weight', sparse = F)
mat[mat == 0] = 1
mat[mat == 2] = 0

a = hclust(as.dist(mat), method = "average", members = NULL)
dend = as.dendrogram(a)
col_aa_red <- ifelse(grepl("Hviridissima", labels(dend)),
                     brewer.pal(3, 'Dark2')[2],
                     ifelse(grepl("Holigactis", labels(dend)),
                            brewer.pal(3, 'Dark2')[1],
                            brewer.pal(3, 'Dark2')[3]))

dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col")

dend2 %>% set("labels_cex", 20) %>% plot(ylab = 'Distance')
#legend(x = 0.7, y = 0.7, legend = c('Hydra oligactis', 'Hydra viridissima', 'Hydra vulgaris'),
#      bty = "n", fill=brewer.pal(3, 'Dark2'), border=NA, cex = 1)
title(main = 'Hierarchical clustering of 1% CR1 repeats', cex.main = 2)
dev.off()

