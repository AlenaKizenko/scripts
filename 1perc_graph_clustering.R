library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(circlize)
library(factoextra)


#-------------------------CR1 1 percent----------------------------------------------------------------------------------

df_hydras_cr1 = read.table('/home/alena/Documents/PhD/CR1_graphs_all_hydras/h_oli_vir_vulchr_blast_1perc_dist_filtered_nonzero.tsv', header = T, sep = '\t')
#df_hydras_10 = read.table('/home/alena/Documents/PhD/CR1_graphs_all_hydras/h_oli_vir_vulchr_blast_1perc_dist_filtered_nonzero.tsv', header = T, sep = '\t')
#grep('Hviridissima', df_hydras_10$Source)
#df_hydras_10 = df_hydras_10[c(1:500, 467800:468300, 146728:146856),]
source = unique(df_hydras_cr1$Source)
target = unique(df_hydras_cr1$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras_cr1 = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hydras_cr1[,1:2], vertices = nodes_hydras_cr1, directed=F) %>%
  set_edge_attr("weight", value = df_hydras_cr1$Weight)

nodes_hydras_cr1$Species = sub("\\_.*", "", nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Nodes = as.character(nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Species = as.factor(nodes_hydras_cr1$Species)
nodes_hydras_cr1$Species = revalue(nodes_hydras_cr1$Species, c(Hvulgarischr = "Hydra vulgaris",
                                                               Holigactis = "Hydra oligactis",
                                                               Hviridissima = 'Hydra viridissima'))

#---------------colors --------------------------------
pal = brewer.pal(3, 'Dark2')
vertex.col = pal[nodes_hydras_cr1$Species]

plot(net,edge.length = E(net)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.color= adjustcolor(vertex.col,alpha = 0.8), vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 3,
     layout=layout_with_fr(net, grid="nogrid")
)
legend('topleft', legend = levels(nodes_hydras_cr1$Species),
       bty = "n", fill=pal, border=NA)


#---------------------clustering----------------------------------------------------
dg <- decompose.graph(net) # returns a list of subgraphs

#-------------------H oligactis subgraph-----------------------------------------
holig = dg[[1]]
df_names = data.frame(Nodes = attr(E(holig), "vnames"))
df_names$Source = sapply(df_names$Nodes, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[[1]]})
df_names$Target = sapply(df_names$Nodes, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[[2]]})
source = unique(df_names$Source)
target = unique(df_names$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras_cr1 = data.frame('Nodes' = c(unique_val, diff1, diff2))

nodes_hydras_cr1$Species = sub("\\_.*", "", nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Nodes = as.character(nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Species = as.factor(nodes_hydras_cr1$Species)
nodes_hydras_cr1$Species = revalue(nodes_hydras_cr1$Species, c(Hvulgarischr = "Hydra vulgaris",
                                                               Holigactis = "Hydra oligactis",
                                                               Hviridissima = 'Hydra viridissima'))

pal = brewer.pal(3, 'Dark2')
vertex.col = pal[nodes_hydras_cr1$Species]

plot(holig,edge.length = E(holig)$weight, edge.arrow.size=.2, edge.curved=0.2,
     vertex.color= adjustcolor(vertex.col,alpha = 0.8), vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 3,
     layout=layout_with_fr(holig, grid="nogrid")
)
legend('topleft', legend = levels(nodes_hydras_cr1$Species),
       bty = "n", fill=pal, border=NA)

write(unique_val, '/home/alena/Documents/PhD/CR1_graphs_all_hydras/Holigactis_subgraph.txt')
#-------------------H vulgaris subgraph-----------------------------------------
hvulg = dg[[61]]

df_names = data.frame(Nodes = attr(E(hvulg), "vnames"))
df_names$Source = sapply(df_names$Nodes, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[[1]]})
df_names$Target = sapply(df_names$Nodes, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[[2]]})
source = unique(df_names$Source)
target = unique(df_names$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras_cr1 = data.frame('Nodes' = c(unique_val, diff1, diff2))

nodes_hydras_cr1$Species = sub("\\_.*", "", nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Nodes = as.character(nodes_hydras_cr1$Nodes)
nodes_hydras_cr1$Species = as.factor(nodes_hydras_cr1$Species)
nodes_hydras_cr1$Species = revalue(nodes_hydras_cr1$Species, c(Hvulgarischr = "Hydra vulgaris",
                                                               Holigactis = "Hydra oligactis",
                                                               Hviridissima = 'Hydra viridissima'))

pal = brewer.pal(3, 'Dark2')
vertex.col = pal[nodes_hydras_cr1$Species]

plot(hvulg,edge.length = E(hvulg)$weight, edge.arrow.size=.2, edge.curved=0.2,
     vertex.color= adjustcolor(vertex.col,alpha = 0.8), vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 3,
     layout=layout_with_fr(hvulg, grid="nogrid")
)
legend('topleft', legend = levels(nodes_hydras_cr1$Species),
       bty = "n", fill=pal, border=NA)


#--------------------rest repeats--------------------------------------

vector = c()

rest = dg[c(2:60, 62:179)]

for (i in 1:177) {
  vector = c(vector, attr(E(rest[[i]]), "vnames"))
}


vector1 = sapply(vector, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[1]})
vector2 = sapply(vector, FUN = function(x) {unlist(strsplit(as.character(x),'\\|'))[2]})

source = unique(vector1)
target = unique(vector2)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
unique_val = unique_val[grep('Holigactis', unique_val)]

write(unique_val, '/home/alena/Documents/PhD/CR1_graphs_all_hydras/Holigactis_other.txt')
