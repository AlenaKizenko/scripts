library(igraph)
library(plyr)
library(RColorBrewer)


#----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
df_hydras_5 = read.table('/Users/alenakizenko/Documents/PhD/project/CR1_graph_fr/hvir_holi_hvul_CR1_fr_blast_out_dist_filtered.tsv',
                            header = T, sep = '\t', stringsAsFactors = F)
df_hydras_5 = subset(df_hydras_5, Source != Target)
source = unique(df_hydras_5$Source)
target = unique(df_hydras_5$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hydras_5[,1:2], vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hydras_5$Weight)
nodes_hydras$Species = sub("\\_.*", "", nodes_hydras$Nodes)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
nodes_hydras$Species = as.factor(nodes_hydras$Species)
nodes_hydras$Species = revalue(nodes_hydras$Species, c(hvir = "Hydra viridissima",
                                                       hvul = "Hydra vulgaris",
                                                       holi = "Hydra oligactis"))
pal = brewer.pal(3, 'Dark2')
vertex.col = pal[nodes_hydras$Species]
#--------------------------subgraph extraction---------------------------------------------
  
clstr_size = data.frame(cluster_size = clusters(net)$csize, clstr_num = c(1:length(clusters(net)$csize)))
clstr_size <- clstr_size[order(-clstr_size$cluster_size),]
  
first = clstr_size[1,2]

membership = as.data.frame(clusters(net)$membership)
membership$Name = rownames(membership)
colnames(membership) = c('Cluster', 'Name')
new = subset(membership, Cluster == first)[,2]
old = subset(membership, Cluster != first)[,2]

write(new, "/Users/alenakizenko/Documents/PhD/CR1_graph_fr/new_repeats.txt")
write(old, "/Users/alenakizenko/Documents/PhD/CR1_graph_fr/old_repeats.txt")

#--------decomposed graph plot-----------
graphs <- decompose.graph(net)
largest <- which.max(sapply(graphs, vcount))

pdf("/Users/alenakizenko/Documents/PhD/project/pictures/rplot.pdf")
plot(net,edge.length = E(net)$weight,
     edge.arrow.size=.2,
     edge.curved=0.2,
     vertex.color = ifelse(V(net)$name %in% new, 'red', 'white'),
     #vertex.color= adjustcolor(vertex.col,alpha = 0.8),
     vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 3,
     layout=layout_with_fr(net, grid="nogrid")
)
dev.off()

#----------tree RVT1----------------------
library(ggtree)
library(RColorBrewer)
library(ape)
library(ggplot2)

#----------------------------CR1 tree 5% of H. oligactis, H. vulgaris, H. viridissima library-------------------------------------------
tree <- read.tree('/Users/alenakizenko/Documents/PhD/CR1_graph_fr/hvir_holi_hvul_5perc_rvt1_align.fasta_prefix.fasta')
tip_labels = tree$tip.label
df = as.data.frame(tip_labels)
tree$Nnode

#tree = root(tree, 'old_hvul_LINE/CR1_rnd-1_family-442_531')
#grp <- list(
 # 'Hydra vulgaris' = tip_labels[grepl("hvul", tip_labels)],
#  'Old' = tip_labels[grepl("old", tip_labels)],
 # 'Old' = tip_labels[grepl("hvir", tip_labels)],
#  'Old' = tip_labels[grepl("holi_LINE/CR1_rnd-1_family-61_586", tip_labels)])


#tree2 = groupOTU(tree, grp, group_name = 'expansion')

grp2 <- list(
  'Hydra vulgaris (new repeats)' = tip_labels[grepl("new_hvul", tip_labels)],
  'Hydra vulgaris (old repeats)' = tip_labels[grepl("old_hvul", tip_labels)],
  'Hydra oligactis (new repeats)' = tip_labels[grepl("new_holi", tip_labels)],
  'Hydra oligactis (old repeats)' = tip_labels[grepl("old_holi", tip_labels)],
  'Hydra oligactis (old repeats)' = tip_labels[grepl("holi_LINE/CR1_rnd-1_family-61_586", tip_labels)],
  'Hydra viridissima (old repeats)' = tip_labels[grepl("hvir", tip_labels)])

tree3 = groupOTU(tree, grp2, group_name = 'species')

p = ggtree(tree3,
           layout = 'circular',
           color = 'gray7'
) +
  # geom_text(aes(label=node), hjust=-.3) +
  geom_point2(aes(color=species, shape = species), size=2.5, alpha = .8, show.legend = T) +
  scale_colour_manual(values=c(brewer.pal(11, 'PuOr')[1], brewer.pal(11, 'PuOr')[4],
                               brewer.pal(11, 'PRGn')[8],
                                          brewer.pal(11, 'PuOr')[11], brewer.pal(11, 'PuOr')[8]
                                          ),
                     # labels = c('Hydra oligactis (new repeats)', 'Hydra oligactis (old repeats)',
                      #           'Hydra vulgaris (old repeats)', 'Hydra vulgaris (new repeats)',
                       #          'Hydra viridissima (old repeats)'),
                      name = '') +
  scale_shape_manual(values = c(16, 17, 17, 16, 17), name = '') +
  theme(
    plot.title = element_blank(),
    # plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/RT_cds_hvir_hvulg_tree_circ.tiff',
       width = 10, height = 10, device = 'tiff')

#----------------------------CR1 tree 5% of H. vulgaris and H. viridissima library-------------------------------------------
tree <- read.tree('/Users/alenakizenko/Documents/PhD/project/rvt_tree/hvir_hvul_rvt1_5perc_cds_align.fasta_prefix.fasta')
tip_labels = tree$tip.label
df = as.data.frame(tip_labels)
tree$tip.label[grepl('hvir', tree$tip.label)]

tree = root(tree, 'hvir_LINE/CR1_rnd-1_family-446_9')

grp2 <- list(
  'Hydra vulgaris (new repeats)' = tip_labels[grepl("new", tip_labels)],
  'Hydra vulgaris (old repeats)' = tip_labels[grepl("old", tip_labels)],
  'Hydra viridissima' = tip_labels[grepl("hvir", tip_labels)])

tree3 = groupOTU(tree, grp2, group_name = 'species')

p = ggtree(tree3,
           layout = 'circular',
           color = 'gray7'
) +
  # geom_text(aes(label=node), hjust=-.3) +
  geom_point2(aes(color=species, shape = species), size=2.5, alpha = .8, show.legend = T) +
  scale_colour_manual(values=c(brewer.pal(11, 'PRGn')[8],
                               brewer.pal(11, 'PuOr')[11], brewer.pal(11, 'PuOr')[8]
  ),
  # labels = c('Hydra oligactis (new repeats)', 'Hydra oligactis (old repeats)',
  #           'Hydra vulgaris (old repeats)', 'Hydra vulgaris (new repeats)',
  #          'Hydra viridissima (old repeats)'),
  name = '') +
  scale_shape_manual(values = c(16, 17, 17, 16, 17), name = '') +
  theme(
    plot.title = element_blank(),
    # plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_cds_hvir_hvul_tree_circ.tiff',
       width = 10, height = 10, device = 'tiff')

