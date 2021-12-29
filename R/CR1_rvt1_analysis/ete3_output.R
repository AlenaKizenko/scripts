library(ggplot2)
library(ggtree)

#-------hydra vulgaris-----------
ete3_df_hvul = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/b_free_b_neut.txt',
                      header = F)
b_neut_hvul = subset(ete3_df_hvul, V1 == 'b_neut')
b_free_hvul = subset(ete3_df_hvul, V1 == 'b_free')

df_hvul = merge(b_free_hvul, b_neut_hvul, by = 'V2')
str(df_hvul)
df_hvul$V3.x = as.numeric(df_hvul$V3.x)
df_hvul$V3.y = as.numeric(df_hvul$V3.y)
df_hvul$lr = 2*(df_hvul$V3.x-df_hvul$V3.y)

df_hvul$pval = 1 - pchisq(df_hvul$lr, df_hvul$V4.x - df_hvul$V4.y)
df_hvul$Selection = sapply(df_hvul$V5.x, function(x) {ifelse(x > 1, 'positive', 'negative')})
df_hvul$branch = sapply(df_hvul$V2, function(x){strsplit(x, '-')[[1]][1]})
df_hvul$branch = as.numeric(df_hvul$branch)

plot_hvul = ggplot(df_hvul, aes(branch, -log10(pval), fill = Selection)) +
  geom_bar(stat = 'identity', color = 'gray7', alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = c(0.87, 0.8),
        legend.text = element_text(size = 18, color = 'black'),
        legend.title = element_text(size = 18, color = 'black'),
        legend.background = element_rect('transparent'),
        panel.grid.major = element_blank()
  ) +
  labs(x = 'Branches', y = '-log10(p-value)') +
  scale_y_continuous(breaks=seq(0, 12, 2), limits=c(0, 14)) +
  scale_x_continuous(breaks=seq(1, 106, 1), limits=c(0, 107)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2]))

ggsave(plot = plot_hvul, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvul_bfreeneut.tiff',
       device = 'tiff', width = 10, height = 10)
  

#-------hydra viridissima-----------
ete3_df_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_viridissima/rvt1_analysis/ete3_b_free_b_neut/result/ete3_b_free_b_neut.txt',
                      header = F)
b_neut_hvir = subset(ete3_df_hvir, V1 == 'b_neut')
b_free_hvir = subset(ete3_df_hvir, V1 == 'b_free')

df_hvir = merge(b_free_hvir, b_neut_hvir, by = 'V2')
str(df_hvir)
df_hvir$V3.x = as.numeric(df_hvir$V3.x)
df_hvir$V3.y = as.numeric(df_hvir$V3.y)
df_hvir$lr = 2*(df_hvir$V3.x-df_hvir$V3.y)

df_hvir$pval = 1 - pchisq(df_hvir$lr, df_hvir$V4.x - df_hvir$V4.y)
df_hvir$Selection = sapply(df_hvir$V5.x, function(x) {ifelse(x > 1, 'positive', 'negative')})
df_hvir$branch = sapply(df_hvir$V2, function(x){strsplit(x, '-')[[1]][1]})
df_hvir$branch = as.numeric(df_hvir$branch)

plot_hvir = ggplot(df_hvir, aes(branch, -log10(pval), fill = Selection)) +
  geom_bar(stat = 'identity', color = 'gray7', alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    #axis.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=18,color = "black"),
    axis.title.x = element_text(size=18,color = "black"),
    axis.title.y = element_text(size=18,color = "black"),
    #strip.text.x = element_text(size = 18,colour = "black"),
    legend.position = c(0.87, 0.8),
    legend.text = element_text(size = 18, color = 'black'),
    legend.title = element_text(size = 18, color = 'black'),
    legend.background = element_rect('transparent')
  ) +
  labs(x = 'Branches', y = '-log10(p-value)') +
  scale_y_continuous(breaks=seq(0, 12, 2), limits=c(0, 14)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2]))

ggsave(plot = plot_hvir, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvir_bfreeneut.tiff',
       device = 'tiff', width = 10, height = 10)

