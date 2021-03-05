library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)

#-------------------------expansion------------------------------------------------------------

genesOfInterest = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/Penelope_repeat_genes_selectedcol.tsv',
                           header = F, sep = ' ')

hist_expansion = ggplot(genesOfInterest, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest <- as.character(unique(genesOfInterest$V4))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


#resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")


test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(myGOdata, test.stat)
pvalFis = score(resultFisher)
hist(pvalFis, xlab = "p-values")

allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
str(allRes)
allRes$classicFisher = as.numeric(allRes$classicFisher)

plot_go = ggplot(allRes, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 18, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO of nearest to brown hydras\' repeat expansion genes',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes$Significant))


ggsave(plot = plot_go, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation/expansion_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')