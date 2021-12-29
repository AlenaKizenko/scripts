library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(gg.gap)
library(ggtree)
library(ape)

bfree_tree_hvir = read.tree('/Users/alenakizenko/Documents/PhD/project/green_hydra_repeat_ann/rvt1_analysis/ete3_b_free_b_neut/output_tree_sign')

tip_labels = bfree_tree_hvir$tip.label
#fb_tree_hvir = root(fb_tree_hvir, 'hvir_LINE/CR1_rnd-1_family-457_56')

#fb_tree_hvir$edge.length = ifelse(fb_tree_hvir$edge.length <= 1, 'positive', 'negative')

p = ggtree(bfree_tree_hvir,
           #  layout = 'circular',
           aes(color=log2(branch.length)),
           branch.length='none'
) +
  geom_tippoint(color = brewer.pal(3, 'Dark2')[1],size=3, alpha = .8) +
  geom_text(aes(label=branch.length), size = 3, color = 'black') +
  scale_color_gradient2(low = brewer.pal(9, name = "Set1")[2],
                        mid = 'white',
                        high = brewer.pal(9, name = "Set1")[1],
                        midpoint =-2,
  ) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p

msaplot(p=p,
        fasta="/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/hvir_rvt1_CR1_codon_align.fasta",
        window=c(200, 300),
        color = c('gray88', brewer.pal(8, 'Accent')[1:4]))


