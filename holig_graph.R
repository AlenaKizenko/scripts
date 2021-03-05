library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(pheatmap)


#-------------------------------green hydra repeats as graph--------------------------------------------------------------------------------------------

df_holig = read.table('/home/alena/Documents/PhD/hydra_oligactis/graph/h_oli_blast_out_dist_notzero.tsv', header = T, sep = '\t')
df_holig$Source = as.factor(df_holig$Source)
df_holig$Target = as.factor(df_holig$Target)
#df_holig = df_holig[1:1000,]
#df_holig = droplevels(df_holig)
nodes = data.frame('id' = levels(unlist(df_holig[,1:2])))
edges = df_holig[, 1:2]
colnames(edges) = c('from', 'to')
edges$length = c(df_holig$Weight)

graph = graph_from_data_frame(d=edges, vertices = nodes, directed = F) %>%
  set_edge_attr("weight", value = edges$length)

pdf(file = '/scratch/kizenko_hydra_project/hydra_oligactis/merged_libs/blast_holig/holig_CR1_graph.pdf')
plot(graph, edge.length = E(graph)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.frame.color="#555555", vertex.color = brewer.pal(3, 'Dark2')[1],
     vertex.size=4, edge.color="gray7", vertex.label=NA, layout=layout.fruchterman.reingold)
dev.off()

#-------------------------------green hydra repeats as heatmap------------------------------------------------------------------------------------


v <- levels(unlist(df_holig[,1:2]))
n <- length(v)                        # number of vertices
e <- matrix(match(as.character(unlist(df_holig[,1:2])), v),ncol=2) # edge list
w <- df_holig$Weight
M <- matrix(0, n, n)                  # set up a distance matrix
M[e] <- w                             # fill it in with edge weights
M <- M + t(M)                         # make this symmetric
dimnames(M) <- list(v, v)             # label the vertices

max(M)

M[M==0] = 1
diag(M) = 0

mybreaks = c(0.01, 0.04, 0.07, 0.1, 0.13, 0.16, 0.19, 0.22, 0.25, 0.3)

pdf(file = '/scratch/kizenko_hydra_project/hydra_oligactis/merged_libs/blast_holig/holig_CR1_heatmap.pdf')
pheatmap(M, color = c("white", "#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"), border_color = "white",
         show_rownames = F, show_colnames = F, breaks = mybreaks)
dev.off()
