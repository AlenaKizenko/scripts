library(genefilter)
library(DESeq2)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

counts_agerev = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_oligactis/revertants/agerev_output.txt',
                     header = F, sep = ' ')
rownames(counts_agerev) = counts_agerev[,1]
counts_agerev = counts_agerev[, -c(1,8)]
str(counts_agerev)
colnames(counts_agerev) = c('age1', 'age2', 'age3', 'rev1', 'rev2', 'rev3')
ann = data.frame('ID' = colnames(counts_agerev),
                 'Condition' = c('age', 'age', 'age', 'rev', 'rev', 'rev'))
ann$Condition = as.factor(ann$Condition)

TE_df = DESeqDataSetFromMatrix(countData = counts_agerev,
                            colData = ann,
                            design = ~Condition)

TE = DESeq(TE_df, betaPrior=TRUE)


res <- results(TE)
head(results(TE, tidy=TRUE)) #
summary(res)
res <- res[order(res$padj),]
head(res)

#pca
vsdata <- varianceStabilizingTransformation(TE_df, blind=FALSE)
plotPCA(vsdata, intgroup="Condition")

#heatmap
ntd <- normTransform(TE)
select <- order(rowMeans(counts(TE,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- data.frame(Condition = colData(TE)[,'Condition'])
pheatmap(assay(vsdata)[select,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols = TRUE,
         show_colnames = TRUE,
         cellwidth = 20, cellheight = 10, angle_col = 45,
         color = colorRampPalette(c(brewer.pal(9, name = "Set1")[2],
                                    'white',
                                    brewer.pal(9, name = "Set1")[1]))(50),
         fontsize = 14, fontsize_row = 9, fontsize_col = 14,
         clustering_distance_rows = "manhattan"
)

# volcanoplot
res_df = as.data.frame(res)
res_df = na.omit(res_df)
str(res_df)
res_df$id = rownames(res_df)
ifelse(abs(res_df$log2FoldChange) > 0.2 & res_df$padj < 0.05, res_df$id, '')

volcanoplot = ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj),
                                 color = ifelse(res_df$padj > 0.05, 'nos',
                                               ifelse(abs(res_df$log2FoldChange) > 0.2, ifelse(res_df$log2FoldChange > 0.2, 'up', 'down'), 'nos')),
                                 label = ifelse(abs(res_df$log2FoldChange) > 0.2 & res_df$padj < 0.05, res_df$id, '')
                                 )
                     ) +
  geom_point(alpha = 0.8) +
  geom_text(size = 4, vjust = 1.5, hjust = 0.5) +
  geom_vline(xintercept = -0.2, color = 'gray55', linetype = 'dashed') +
  geom_vline(xintercept = 0.2, color = 'gray55', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), color = 'gray55', linetype = 'dashed') +
  xlab("log2FoldChange") +
  ylab("-log10(p-value adjusted)") +
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), limits = c(-0.6, 0.6)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.y = element_text(size=16, color = "black"),
        axis.title.x = element_text(size=16, color = "black")
  ) +
  scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[10], brewer.pal(8, 'Dark2')[8], brewer.pal(11, 'RdYlBu')[1]),
                     name = 'Differential expression')

volcanoplot




