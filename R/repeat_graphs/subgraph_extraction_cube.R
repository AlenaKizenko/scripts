library(igraph)
library(plyr)


#----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
df_hydras = read.table('/scratch/kizenko_hydra_project/blast_all_CR1_fr/hvir_holi_hvul_CR1_fr_blast_out_dist_filt.tsv',
                         header = T, sep = '\t', stringsAsFactors = F)
df_hydras = subset(df_hydras, Source != Target)
source = unique(df_hydras$Source)
target = unique(df_hydras$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hydras[,1:2], vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hydras$Weight)

#--------------------------subgraph extraction---------------------------------------------

clstr_size = data.frame(cluster_size = clusters(net)$csize, clstr_num = c(1:length(clusters(net)$csize)))
clstr_size <- clstr_size[order(-clstr_size$cluster_size),]

first = clstr_size[1,2]

membership = as.data.frame(clusters(net)$membership)
membership$Name = rownames(membership)
colnames(membership) = c('Cluster', 'Name')
new = subset(membership, Cluster == first)[,2]
old = subset(membership, Cluster != first)[,2]

write(new, "/scratch/kizenko_hydra_project/blast_all_CR1_fr/new_repeats_fr.txt")
write(old, "/scratch/kizenko_hydra_project/blast_all_CR1_fr/old_repeats_fr.txt")