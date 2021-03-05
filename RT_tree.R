library(ggtree)
#library(data.tree)
library(RColorBrewer)
library(ape)
#library(phylobase)
library(ggplot2)


#----------------------------CR1 tree 10% of H. oligactis, H. vulgaris, H. viridissima library-------------------------------------------
tree <- read.tree('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/LINE_CR1_repeatsonly_all_prefix_filtered_tree')
tip_labels = tree$tip.label
tree$Nnode

tree = root(tree, 'rest_LINE/CR1_rnd-6_family-545_63')
grp <- list(
            'Expanded' = tip_labels[grepl("expansion", tip_labels)],
            'Non-expanded' = tip_labels[grepl("rest", tip_labels)])


tree2 = groupOTU(tree, grp)
tip_labels[grepl('rest', tip_labels)]

p = ggtree(tree2,
           layout = 'circular'
         #  aes(color=group)
         #  branch.length='none'
         ) +
#  geom_text(aes(label=node), hjust=-.3) +
  geom_tippoint(aes(color=group, shape = group), size=1.5, alpha = .8) +
  scale_colour_manual(name = 'Repeat type', values=c(brewer.pal(4, 'Set2')[3], brewer.pal(4, 'Set2')[4])) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face="bold"),
  # plot.subtitle = element_text(size = 16, hjust = 0),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(title = 'Phylogeny of reverse transcriptases in CR1 repeats')
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/RT_hvulg_tree.tiff',
       width = 10, height = 10, device = 'tiff')


treeio::tree_subset(p, 1:50)

p2 <- p %>% collapse(node=21) + 
  geom_point2(aes(subset=(node==21)), shape=21, size=5, fill='green')
p2 <- collapse(p2, node=23) + 
  geom_point2(aes(subset=(node==23)), shape=23, size=5, fill='red')
print(p2)
expand(p2, node=23) %>% expand(node=21)