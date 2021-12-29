library(ggplot2)
library(RColorBrewer)


#------------read GC content--------------------------------
sym = read.csv2('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/aposym/1695c_sym1_1_paired_counts.tab',
                sep = '\t', header = F)

apo = read.csv2('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/aposym/1695c_apo1_1_paired_counts.tab',
                sep = '\t', header = F)

sym$V1 = as.numeric(sym$V1)
apo$V1 = as.numeric(apo$V1)
sym$name = rep('Symbiotic',each=length(sym))
colnames(sym) = c('Read GC content', 'Hydra')
apo$name = rep('Aposymbiotic',each=length(apo))
colnames(apo) = c('Read GC content', 'Hydra')

df3 = rbind(sym[1:1000,], apo[1:1000,])

sym_plot = ggplot(df3, aes(Hydra, `Read GC content`, fill = Hydra)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, size = 1) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=18),
    axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
    axis.text.y = element_text(size = 18),
    axis.title=element_text(size=18),
    legend.position = 'none'
  ) +
  scale_fill_manual(values = c(brewer.pal(11, 'PRGn')[3],brewer.pal(11, 'PRGn')[9]))

ggsave(plot = sym_plot,
        filename = '/Users/alenakizenko/Documents/PhD/pictures/plot_GC_apo_sym.tiff',
        width = 10, height = 10, device = 'tiff')

#------------gene DE analysis---------------------
library(DESeq2)


df_TE = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_viridissima/aposym/apo_sym_counts.txt',
                    skip = 1)

rownames(df_TE) = df_TE[,1]
df_TE = df_TE[, -c(1,2,3,4,5,6)]
colnames(df_TE) = c("1695c_apo", "1695c_sym", "CA07a_apo",
                    "CA07a_sym", "GBR01b_apo", "GBR01b_sym",
                    "PAN04a_apo", "PAN04a_sym")

ann = data.frame('ID' = colnames(df_TE),
                 'Condition' = c('apo', 'sym',
                                 'apo', 'sym',
                                 'apo', 'sym',
                                 'apo', 'sym'),
                 'Strain' = c("1695c", "1695c",
                              "CA07a", "CA07a",
                              "GBR01b", "GBR01b",
                              "PAN04a", "PAN04a")
)
ann$ID = as.factor(ann$ID)
ann$Condition = as.factor(ann$Condition)
ann$Strain = as.factor(ann$Strain)

df_deseq2 = DESeqDataSetFromMatrix(countData =df_TE, colData = ann, design = ~Strain)
condition_res = DESeq(df_deseq2)
res <- results(condition_res)
head(results(condition_res, tidy=TRUE)) #
summary(res)
res <- res[order(res$padj),]
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

vsdata <- vst(df_deseq2, blind=FALSE)
plotPCA(vsdata, intgroup="Strain")

