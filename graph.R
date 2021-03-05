library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)

#----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
df_hydras_10 = read.table('/home/alena/Documents/PhD/CR1_graphs_all_hydras/h_olig_vir_vulg_distances_filtered.tsv',
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

#pdf(file = '/scratch/kizenko_hydra_project/three_hydras_blastn_10perc/hydras10perc_CR1.pdf')
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
title(main = 'H. viridissima, H. oligactis, H. vulgaris LINE/CR1 repeat graph', cex.main = 1.4)
dev.off()


#------------------------------------------Hydra viridissima---------------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(network)
library(pheatmap)


df_hvir = read.table('/home/alena/Documents/PhD/green_hydra_repeat_ann/graph_approach/h_vir_blast_out_dist_filt.tsv',
                     header = T, sep = '\t')
df_hvir = subset(df_hvir, Source != Target)
source = unique(df_hvir$Source)
target = unique(df_hvir$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))

net <- graph_from_data_frame(d=df_hvir, vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hvir$Weight)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)

plot(net,edge.length = E(net)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.frame.color= "#555555",
     vertex.label=NA, vertex.color= adjustcolor(brewer.pal(3, 'Dark2')[1],alpha = 0.7),
     edge.color="gray8", vertex.size = 2,
     layout=layout_with_fr(net, grid="nogrid"))
title(main = 'Hydra viridissima CR1 repeat graph', cex.main = 2)
dev.off()

#----------------------------------------hydra vulgaris-----------------------------------------------------------------------------------------------

library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)


df_hvulg = read.table('/scratch/kizenko_hydra_project/hydra_vulgaris_chrom/blast_CR1_hvulgchr/h_vulgchr_blast_out_dist_filt.tsv', header = T, sep = '\t')
df_hvulg = subset(df_hvulg, Source != Target)
source = unique(df_hvulg$Source)
target = unique(df_hvulg$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))
net <- graph_from_data_frame(d=df_hvulg, vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_hvulg$Weight)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
pdf(file = '/scratch/kizenko_hydra_project/hydra_vulgaris_chrom/blast_CR1_hvulgchr/hvulgchr_CR1_graph.pdf')
plot(net,edge.length = E(net)$weight, edge.arrow.size=.2, edge.curved=0.2, vertex.frame.color= "#555555",
     vertex.label=NA, vertex.color= adjustcolor(brewer.pal(3, 'Dark2')[3],alpha = 0.7),
     edge.color="gray8", vertex.size = 2,
     layout=layout_with_fr(net, grid="nogrid"))
title(main = 'Hydra vulgaris CR1 repeat graph', cex.main = 2)
dev.off()


#--------------------------------hydra oligactis--------------------------------------------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)


df_holig = read.table('/scratch/kizenko_hydra_project/hydra_oligactis/merged_libs/blast_holig/h_oli_blast_out_dist_filt.tsv', header = T, sep = '\t')
df_holig = subset(df_holig, Source != Target)
source = unique(df_holig$Source)
target = unique(df_holig$Target)
diff1 = setdiff(source, target)
diff2 = setdiff(target, source)
unique_val = as.vector(intersect(source, target))
nodes_hydras = data.frame('Nodes' = c(unique_val, diff1, diff2))

net <- graph_from_data_frame(d=df_holig, vertices = nodes_hydras, directed=F) %>%
  set_edge_attr("weight", value = df_holig$Weight)

nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
pdf(file = '/scratch/kizenko_hydra_project/hydra_oligactis/merged_libs/blast_holig/holig_CR1_graph.pdf')
plot(net,edge.length = E(net)$weight, edge.arrow.size=.2, edge.curved=0.2,
     vertex.color= adjustcolor(brewer.pal(3, 'Dark2')[2],alpha = 0.7), vertex.frame.color= "#555555",
     vertex.label=NA,
     edge.color="gray8", vertex.size = 2,
     layout=layout_with_fr(net, grid="nogrid"))
title(main = 'Hydra oligactis CR1 repeat graph', cex.main = 2)
dev.off()


#-------------------------Penelope 1 percent----------------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)

df_hydras_penelope1 = read.table('/scratch/kizenko_hydra_project/blast_all_Penelope/h_oli_vir_vul_Penelope_distances_filtered.tsv',
                                 header = T, sep = '\t')

df_hydras_penelope1 = subset(df_hydras_penelope1, Source != Target)

cnt_hvulg_in = 0

cnt_hvulg_out = 0

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Hvulgaris', df_hydras_penelope1[row, 1]) == grepl('Hvulgaris', df_hydras_penelope1[row, 2])
     & grepl('Hvulgaris', df_hydras_penelope1[row, 1]) == TRUE){
    cnt_hvulg_in = cnt_hvulg_in + 1
  }
}

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Hvulgaris', df_hydras_penelope1[row, 1]) != grepl('Hvulgaris', df_hydras_penelope1[row, 2])
     & (grepl('Hvulgaris', df_hydras_penelope1[row, 1]) == TRUE
        | grepl('Hvulgaris', df_hydras_penelope1[row, 2])) == TRUE){
    cnt_hvulg_out = cnt_hvulg_out + 1
  }
}
# count all Hydra oligactis inner and outer repeats

cnt_holig_in = 0

cnt_holig_out = 0

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Holigactis', df_hydras_penelope1[row, 1]) == grepl('Holigactis', df_hydras_penelope1[row, 2])
     & grepl('Holigactis', df_hydras_penelope1[row, 1]) == TRUE){
    cnt_holig_in = cnt_holig_in + 1
  }
}

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Holigactis', df_hydras_penelope1[row, 1]) != grepl('Holigactis', df_hydras_penelope1[row, 2])
     & (grepl('Holigactis', df_hydras_penelope1[row, 1]) == TRUE
        | grepl('Holigactis', df_hydras_penelope1[row, 2])) == TRUE){
    cnt_holig_out = cnt_holig_out + 1
  }
}


# count all Hydra viridissima inner and outer repeats

cnt_hvir_in = 0

cnt_hvir_out = 0

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Hviridissima', df_hydras_penelope1[row, 1]) == grepl('Hviridissima', df_hydras_penelope1[row, 2])
     & grepl('Hviridissima', df_hydras_penelope1[row, 1]) == TRUE){
    cnt_hvir_in = cnt_hvir_in + 1
  }
}

for (row in 1:nrow(df_hydras_penelope1)){
  if(grepl('Hviridissima', df_hydras_penelope1[row, 1]) != grepl('Hviridissima', df_hydras_penelope1[row, 2])
     & (grepl('Hviridissima', df_hydras_penelope1[row, 1]) == TRUE
        | grepl('Hviridissima', df_hydras_penelope1[row, 2])) == TRUE){
    cnt_hvir_out = cnt_hvir_out + 1
  }
}

# barplot Penelope edges

cnt_hvir_in_perc = (cnt_hvir_in *100)/(cnt_hvir_in + cnt_hvir_out)
cnt_hvir_out_perc = (cnt_hvir_out *100)/(cnt_hvir_in + cnt_hvir_out)

cnt_holig_in_perc = (cnt_holig_in *100)/(cnt_holig_in + cnt_holig_out)
cnt_holig_out_perc = (cnt_holig_out *100)/(cnt_holig_in + cnt_holig_out)

cnt_hvulg_in_perc = (cnt_hvulg_in *100)/(cnt_hvulg_in + cnt_hvulg_out)
cnt_hvulg_out_perc = (cnt_hvulg_out *100)/(cnt_hvulg_in + cnt_hvulg_out)


penelope_edges = data.frame('Count' = c(cnt_hvir_in_perc, cnt_hvir_out_perc,
                                   cnt_holig_in_perc, cnt_holig_out_perc,
                                   cnt_hvulg_in_perc, cnt_hvulg_out_perc),
                       'Species' = c('Hydra viridissima', 'Hydra viridissima',
                                     'Hydra oligactis', 'Hydra oligactis',
                                     'Hydra vulgaris', 'Hydra vulgaris'),
                       'Edge type' = c('in', 'out', 'in', 'out', 'in', 'out'))

plot = ggplot(penelope_edges, aes(Edge.type, Count, fill = Edge.type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(~Species) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  labs(y = 'Percentage of edges') +
  scale_fill_manual(values = brewer.pal(11, 'RdYlBu')[c(4,9)],
                    name="Edge type",
                    labels=c("inner", "outer")) +
  scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))
ggsave(plot = plot,
       filename = '/scratch/kizenko_hydra_project/blast_all_Penelope/barplot_Penelope.tiff',
       width = 10, height = 10, device = 'tiff')
#-------------------------CR1 1 percent----------------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)

df_hydras_cr1 = read.table('/scratch/kizenko_hydra_project/blast_all_CR1/h_olig_vir_vulg_CR1_distances_filtered.tsv',
                           header = T, sep = '\t', stringsAsFactors = F)

df_hydras_cr1 = subset(df_hydras_cr1, Source != Target)

cnt_hvulg_in = 0
cnt_hvulg_out = 0

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Hvulgaris', df_hydras_cr1[row, 1]) == grepl('Hvulgaris', df_hydras_cr1[row, 2])
     & grepl('Hvulgaris', df_hydras_cr1[row, 1]) == TRUE){
    cnt_hvulg_in = cnt_hvulg_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Hvulgaris', df_hydras_cr1[row, 1]) != grepl('Hvulgaris', df_hydras_cr1[row, 2])
     & (grepl('Hvulgaris', df_hydras_cr1[row, 1]) == TRUE
        | grepl('Hvulgaris', df_hydras_cr1[row, 2])) == TRUE){
    cnt_hvulg_out = cnt_hvulg_out + 1
  }
}
# count all Hydra oligactis inner and outer repeats

cnt_holig_in = 0

cnt_holig_out = 0

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Holigactis', df_hydras_cr1[row, 1]) == grepl('Holigactis', df_hydras_cr1[row, 2])
     & grepl('Holigactis', df_hydras_cr1[row, 1]) == TRUE){
    cnt_holig_in = cnt_holig_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Holigactis', df_hydras_cr1[row, 1]) != grepl('Holigactis', df_hydras_cr1[row, 2])
     & (grepl('Holigactis', df_hydras_cr1[row, 1]) == TRUE
     | grepl('Holigactis', df_hydras_cr1[row, 2])) == TRUE){
    cnt_holig_out = cnt_holig_out + 1
  }
}


# count all Hydra viridissima inner and outer repeats

cnt_hvir_in = 0

cnt_hvir_out = 0

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Hviridissima', df_hydras_cr1[row, 1]) == grepl('Hviridissima', df_hydras_cr1[row, 2])
     & grepl('Hviridissima', df_hydras_cr1[row, 1]) == TRUE){
    cnt_hvir_in = cnt_hvir_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('Hviridissima', df_hydras_cr1[row, 1]) != grepl('Hviridissima', df_hydras_cr1[row, 2])
     & (grepl('Hviridissima', df_hydras_cr1[row, 1]) == TRUE
    | grepl('Hviridissima', df_hydras_cr1[row, 2])) == TRUE){
    cnt_hvir_out = cnt_hvir_out + 1
  }
}

# barplot CR1 edges

cnt_hvir_in_perc = (cnt_hvir_in *100)/(cnt_hvir_in + cnt_hvir_out)
cnt_hvir_out_perc = (cnt_hvir_out *100)/(cnt_hvir_in + cnt_hvir_out)

cnt_holig_in_perc = (cnt_holig_in *100)/(cnt_holig_in + cnt_holig_out)
cnt_holig_out_perc = (cnt_holig_out *100)/(cnt_holig_in + cnt_holig_out)

cnt_hvulg_in_perc = (cnt_hvulg_in *100)/(cnt_hvulg_in + cnt_hvulg_out)
cnt_hvulg_out_perc = (cnt_hvulg_out *100)/(cnt_hvulg_in + cnt_hvulg_out)


cr1_edges = data.frame('Count' = c(cnt_hvir_in_perc, cnt_hvir_out_perc,
                                   cnt_holig_in_perc, cnt_holig_out_perc,
                                   cnt_hvulg_in_perc, cnt_hvulg_out_perc),
                       'Species' = c('Hydra viridissima', 'Hydra viridissima',
                                     'Hydra oligactis', 'Hydra oligactis',
                                     'Hydra vulgaris', 'Hydra vulgaris'),
                       'Edge type' = c('in', 'out', 'in', 'out', 'in', 'out'))

plot = ggplot(cr1_edges, aes(Edge.type, Count, fill = Edge.type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(~Species) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    strip.text.x = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  labs(y = 'Percentage of edges') +
  scale_fill_manual(values = brewer.pal(11, 'RdYlBu')[c(4,9)],
                    name="Edge type",
                    labels=c("inner", "outer")) +
  scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))
ggsave(plot = plot,
       filename = '/scratch/kizenko_hydra_project/blast_all_CR1/barplot_CR1.tiff',
       width = 10, height = 10, device = 'tiff')


