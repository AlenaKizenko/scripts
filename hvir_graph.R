library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(pheatmap)


#-------------------------------green hydra repeats as graph--------------------------------------------------------------------------------------------

df_hvir = read.table('/home/alena/Documents/PhD/green_hydra_repeat_ann/graph_approach/h_vir_blast_out_dist_filt_notzero.tsv', header = T, sep = '\t')
nodes = data.frame('id' = levels(unlist(df_hvir[,1:2])))
edges = df_hvir[, 1:2]
colnames(edges) = c('from', 'to')
edges$length = c(df_hvir$Weight)

graph = graph_from_data_frame(d=edges[1:2], vertices = nodes, directed = F) %>%
  set_edge_attr("weight", value = edges$length)

plot(graph, edge.length = E(graph)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.frame.color="#555555", vertex.color = brewer.pal(3, 'Dark2')[2],
     vertex.size=4, edge.color="gray7", vertex.label=NA, layout=layout.fruchterman.reingold)
#dev.off()

#-------------------------------green hydra repeats as heatmap------------------------------------------------------------------------------------


v <- levels(unlist(df_hvir[,1:2]))

n <- length(v)                        # number of vertices
e <- matrix(match(as.character(unlist(df_hvir[,1:2])), v),ncol=2) # edge list
w <- df_hvir$Weight
M <- matrix(0, n, n)                  # set up a distance matrix
M[e] <- w                             # fill it in with edge weights
M <- M + t(M)                         # make this symmetric
dimnames(M) <- list(v, v)             # label the vertices

M[M==0] = 1

diag(M) = 0
mybreaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)

heatmap_hvir = pheatmap(M, color = rev(c("white", "#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#8C2D04")), border_color = "white",
         show_rownames = F, show_colnames = F, breaks = mybreaks)

df_clusters = as.data.frame(cutree(heatmap_hvir$tree_row, k = 15))
#df_clusters$`cutree(heatmap_hvir$tree_row, k = 15)` = as.character(df_clusters$`cutree(heatmap_hvir$tree_row, k = 15)`)
str(df_clusters)
ggplot(df_clusters, aes(df_clusters$`cutree(heatmap_hvir$tree_row, k = 15)`)) +
  geom_histogram(stat = 'count', fill = "#FDAE6B") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12, 13, 14, 15)) +
  xlab('Cluster') +
  ylab('Number of repeats') +
  theme_bw() +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18))

