library(igraph)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggplotify)
library(svglite)

#----------------H. viridissima, H. oligactis, H. vulgaris CR1 repeat graph----------------------------------------------------------
df_hydras_5 = read.table('/Users/alenakizenko/Documents/PhD/CR1_graph_fr/hvir_holi_hvul_CR1_fr_blast_out_dist_filtered.tsv',
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
str(nodes_hydras)
nodes_hydras$Nodes = as.character(nodes_hydras$Nodes)
nodes_hydras$Species = as.factor(nodes_hydras$Species)
nodes_hydras$Species
nodes_hydras$Species = revalue(nodes_hydras$Species, c(hvul = "Hydra vulgaris",
                                                       holi = "Hydra oligactis",
                                                       hvir = 'Hydra viridissima'))

# colors
pal =c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1], brewer.pal(3, 'Dark2')[3])
vertex.col = pal[nodes_hydras$Species]

svglite('/Users/alenakizenko/Documents/PhD/pictures/article/hydras5perc_CR1.svg')
plot_cr1 = plot(net,edge.length = E(net)$weight,
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
       bty = "n", fill=pal, border=NA, cex = 1)
invisible(dev.off())

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
#-------------------------CR1----------------------------------------------------------------------------------
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
  if(grepl('hvul', df_hydras_cr1[row, 1]) == grepl('hvul', df_hydras_cr1[row, 2])
     & grepl('hvul', df_hydras_cr1[row, 1]) == TRUE){
    cnt_hvulg_in = cnt_hvulg_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('hvul', df_hydras_cr1[row, 1]) != grepl('hvul', df_hydras_cr1[row, 2])
     & (grepl('hvul', df_hydras_cr1[row, 1]) == TRUE
        | grepl('hvul', df_hydras_cr1[row, 2])) == TRUE){
    cnt_hvulg_out = cnt_hvulg_out + 1
  }
}
# count all Hydra oligactis inner and outer repeats

cnt_holig_in = 0

cnt_holig_out = 0

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('holi', df_hydras_cr1[row, 1]) == grepl('holi', df_hydras_cr1[row, 2])
     & grepl('holi', df_hydras_cr1[row, 1]) == TRUE){
    cnt_holig_in = cnt_holig_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('holi', df_hydras_cr1[row, 1]) != grepl('holi', df_hydras_cr1[row, 2])
     & (grepl('holi', df_hydras_cr1[row, 1]) == TRUE
     | grepl('holi', df_hydras_cr1[row, 2])) == TRUE){
    cnt_holig_out = cnt_holig_out + 1
  }
}


# count all Hydra viridissima inner and outer repeats

cnt_hvir_in = 0

cnt_hvir_out = 0

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('hvir', df_hydras_cr1[row, 1]) == grepl('hvir', df_hydras_cr1[row, 2])
     & grepl('hvir', df_hydras_cr1[row, 1]) == TRUE){
    cnt_hvir_in = cnt_hvir_in + 1
  }
}

for (row in 1:nrow(df_hydras_cr1)){
  if(grepl('hvir', df_hydras_cr1[row, 1]) != grepl('hvir', df_hydras_cr1[row, 2])
     & (grepl('hvir', df_hydras_cr1[row, 1]) == TRUE
    | grepl('hvir', df_hydras_cr1[row, 2])) == TRUE){
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
       filename = '/scratch/kizenko_hydra_project/blast_all_CR1/barplot_CR1.svg',
       width = 10, height = 10, device = 'svg')


