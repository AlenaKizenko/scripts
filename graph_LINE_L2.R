library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)

#----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
df_hydras_1 = read.table('/home/alena/Documents/PhD/blast_L2_1perc/h_oli_vir_vul_L2_1perc_blast_out_distances_filt.tsv',
                          header = T, sep = '\t', stringsAsFactors = F)
df_hydras_1 = subset(df_hydras_1, Source != Target)
source = unique(df_hydras_1$Source)
target = unique(df_hydras_1$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hydras_1[,1:2], vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hydras_1$Weight)

nodes_hydras$Species = sub("\\_.*", "", nodes_hydras$Nodes)
str(nodes_hydras)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
nodes_hydras$Species = as.factor(nodes_hydras$Species)
nodes_hydras$Species = revalue(nodes_hydras$Species, c(Hvulgaris = "Hydra vulgaris",
                                                       Holigactis = "Hydra oligactis",
                                                       Hviridissima = 'Hydra viridissima'))

# colors
pal =c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1], brewer.pal(3, 'Dark2')[3])
vertex.col = pal[nodes_hydras$Species]

#df(file = '/home/alena/Documents/PhD/blast_L2_1perc/hydras1perc_L2.pdf')
plot(net,edge.length = E(net)$weight,
     edge.arrow.size=.2,
     edge.curved=0.2,
     vertex.color= adjustcolor(vertex.col,alpha = 0.8),
     vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8",
     vertex.size = 3,
     layout=layout_with_fr(net, grid="nogrid")
)
legend('topleft', legend = levels(nodes_hydras$Species),
       bty = "n", fill=pal, border=NA, cex = 1.3)
title(main = 'H. viridissima, H. oligactis, H. vulgaris LINE/L2 repeat graph', cex.main = 1.4)
dev.off()