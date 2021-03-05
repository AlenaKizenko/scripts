library(ggtree)
#library(data.tree)
library(RColorBrewer)
#library(ape)
#library(phylobase)
library(ggplot2)


#----------------------------CR1 tree 10% of H. oligactis, H. vulgaris, H. viridissima library-------------------------------------------
tree <- read.tree('/home/alena/Documents/PhD/tree/H_vir_oli_vulg_tree')
tip_labels = tree$tip.label
tree$Nnode

grp <- list(
            'Hydra oligactis'  = tip_labels[grepl("Holigactis", tip_labels)],
            'Hydra vulgaris' = tip_labels[grepl("Hvulgaris", tip_labels)],
            'Hydra viridissima' = tip_labels[grepl("Hviridissima", tip_labels)])


tree2 = groupOTU(tree, grp)


p = ggtree(tree2, layout = 'circular', aes(color=group), branch.length='none') +
  scale_colour_manual(name = 'Species', values=c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1], brewer.pal(3, 'Dark2')[3])) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face="bold"),
    plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
    ) +
  labs(title = 'CR1 repeats diversity (1% of library)')

gzoom(tree2, 1:50, subtree = T)


#----------------------------Penelope tree 10% of H. oligactis, H. vulgaris, H. viridissima library-------------------------------------------
tree_penelope <- read.tree('/home/alena/Documents/PhD/tree/H_vir_oli_vulchr_tree_Penelope1')
tip_labels_penelope = tree_penelope$tip.label

grp_penelope <- list(
            'Hydra oligactis'  = tip_labels_penelope[grepl("Holigactis", tip_labels_penelope)],
            'Hydra vulgaris' = tip_labels_penelope[grepl("Hvulgaris", tip_labels_penelope)],
            'Hydra viridissima' = tip_labels_penelope[grepl("Hviridissima", tip_labels_penelope)])


tree2_penelope = groupOTU(tree_penelope, grp_penelope)


p_penelope = ggtree(tree2_penelope, layout = 'fan', aes(color=group), branch.length='none') +
  scale_colour_manual(name = 'Species', values=c(brewer.pal(3, 'Dark2')[2], brewer.pal(3, 'Dark2')[1], brewer.pal(3, 'Dark2')[3])) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face="bold"),
    plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  labs(title = 'Penelope repeats diversity (1% of library)')

gzoom(tree2_penelope, 1:50, subtree = T)
viewClade(ggtree(tree2_penelope), MRCA(tree2_penelope, tip=c("Holigactis_21527_LINE/Penelope", "Hvulgarischr_5494_LINE/Penelope")))


#----------------------------RVT1 tree of H. oligactis, H. vulgaris, H. viridissima CR1 library-------------------------------------------
tree_rvt1 <- read.tree('/home/alena/Documents/PhD/rvt1_all_hydras_tree/rvt1_all_hydras_tree')
tip_labels_rvt1 = tree_rvt1$tip.label

grp_rvt1 <- list(
                 'Hydra vulgaris' = tip_labels_rvt1[grepl("Hvulgarischr", tip_labels_rvt1)],
                 'Hydra oligactis'  = tip_labels_rvt1[grepl("Holigactis", tip_labels_rvt1)],
                 'Hydra viridissima' = tip_labels_rvt1[grepl("Hviridissima", tip_labels_rvt1)])


tree2_rvt1 = groupOTU(tree_rvt1, grp_rvt1)


p_rvt1 = ggtree(tree2_rvt1, layout = 'circular', aes(color=group), branch.length='none') +
  scale_colour_manual(name = 'Species', values=brewer.pal(3, 'Dark2')) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face="bold"),
    plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  labs(title = 'CR1 transposons RVT1 tree', subtitle = 'Without branch lengths')

subtree_penelope = gzoom(tree2_penelope, 1:10, subtree = T)
viewClade(ggtree(tree2_penelope), MRCA(tree2_penelope, tip=c("Holigactis_21527_LINE/Penelope", "Hvulgarischr_5494_LINE/Penelope")))


p1_rvt1 = ggtree(tree2_rvt1, layout = 'circular', aes(color=group), size=.8) +
  scale_colour_manual(name = 'Species', values=brewer.pal(3, 'Dark2')) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face="bold"),
    plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  labs(title = 'CR1 transposons RVT1 tree', subtitle = 'With branch lengths')


