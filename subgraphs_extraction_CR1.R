library(igraph)
#library(ggplot2)
#library(RColorBrewer)
library(plyr)

  
  #----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
  df_hydras_10 = read.table('/scratch/kizenko_hydra_project/blast_all_Penelope/h_oli_vir_vul_Penelope_distances_filtered.tsv',
                            header = T, sep = '\t', stringsAsFactors = F)
  df_hydras_10 = subset(df_hydras_10, Source != Target)
  source = unique(df_hydras_10$Source)
  target = unique(df_hydras_10$Target)
  diff1 = setdiff(source, target)
  diff2 = setdiff(target, source)
  unique_val = as.vector(intersect(source, target))
  nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
  net <- graph_from_data_frame(d=df_hydras_10[,1:2], vertices = nodes_hydras, directed=F) %>%
    set_edge_attr("weight", value = df_hydras_10$Weight)
  
  
  #nodes_hydras$Species = sub("\\_.*", "", nodes_hydras$Nodes)
  #str(nodes_hydras)
  #nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
  #nodes_hydras$Species = as.factor(nodes_hydras$Species)
  #nodes_hydras$Species = revalue(nodes_hydras$Species, c(Hvulgaris = "Hydra vulgaris",
   #                                                      Holigactis = "Hydra oligactis",
    #                                                     Hviridissima = 'Hydra viridissima'))
  
  # colors
  #pal =c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1], brewer.pal(3, 'Dark2')[3])
  #vertex.col = pal[nodes_hydras$Species]
  
  #pdf(file = '/scratch/kizenko_hydra_project/blast_all_CR1/graph_CR1.pdf')
  #plot(net,edge.length = E(net)$weight,
   #    edge.arrow.size=.2,
    #   edge.curved=0.2,
     #  vertex.color= adjustcolor(vertex.col,alpha = 0.8),
  #     vertex.frame.color= "#555555",
   #    vertex.label=NA,
    #   edge.color="gray8",
     #  vertex.size = 3,
      # layout=layout_with_fr(net, grid="nogrid")
  #)
  #legend('topleft', legend = levels(nodes_hydras$Species),
   #      bty = "n", fill=pal, border=NA, cex = 1.3)
  #title(main = 'H. viridissima, H. oligactis, H. vulgaris LINE/CR1 repeat graph', cex.main = 1.4)
  #dev.off()
  
  #--------------------------subgraph extraction---------------------------------------------
  
  clstr_size = data.frame(cluster_size = clusters(net)$csize, clstr_num = c(1:length(clusters(net)$csize)))
  clstr_size <- clstr_size[order(-clstr_size$cluster_size),]
  
  first = clstr_size[1,2]
  
  membership = as.data.frame(clusters(net)$membership)
  membership$Name = rownames(membership)
  colnames(membership) = c('Cluster', 'Name')
  hvulgaris = subset(membership, Cluster == first)[,2]
  rest = subset(membership, Cluster != first)[,2]
  write(hvulgaris, "/scratch/kizenko_hydra_project/blast_all_Penelope/expansion_repeats.txt")
  write(rest, "/scratch/kizenko_hydra_project/blast_all_Penelope/rest_rest.txt")
    
