library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)


df_hydras_10 = read.table('/scratch/kizenko_hydra_project/three_hydras_blastn_Penelope_10perc/h_oli_vul_vir_blast_Penelope10_dist_filtered_nonzero.tsv', header = T, sep = '\t')
#df_hydras_10 = read.table('/home/alena/Documents/PhD/Penelope_10perc/h_oli_vul_vir_blast_Penelope10_dist_filtered_nonzero.tsv', header = T, sep = '\t')
#df_hydras_10 = read.table('/home/alena/Documents/PhD/CR1_graphs_all_hydras/h_oli_vir_vulchr_blast_1perc_dist_filtered_nonzero.tsv', header = T, sep = '\t')
#grep('Holigactis', df_hydras_10$Source)
#df_hydras_10 = df_hydras_10[c(1:400, 50897:51397, 50301:50701),]
source = unique(df_hydras_10$Source)
target = unique(df_hydras_10$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hydras_10[,1:2], vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hydras_10$Weight)

nodes_hydras$Species = sub("\\_.*", "", nodes_hydras$Nodes)
str(nodes_hydras)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
nodes_hydras$Species = as.factor(nodes_hydras$Species)
nodes_hydras$Species = revalue(nodes_hydras$Species, c(Hvulgarischr = "Hydra vulgaris", Holigactis = "Hydra oligactis", Hviridissima = 'Hydra viridissima'))

#---------------colors --------------------------------
pal = brewer.pal(3, 'Dark2')
vertex.col = pal[nodes_hydras$Species]

#-------------------------------------------------------------------

pdf(file = '/scratch/kizenko_hydra_project/three_hydras_blastn_Penelope_10perc/hydras10perc_Penelope.pdf')
plot(net,edge.length = E(net)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.color= adjustcolor(vertex.col,alpha = 0.8), vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 3,
     layout=layout_with_fr(net, grid="nogrid")
)
legend('topleft', legend = levels(nodes_hydras$Species),
       bty = "n", fill=pal, border=NA)
dev.off()