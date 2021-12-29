library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)

#-------------------------DNA transposons------------------------------------------------------------

genesOfInterest_DNA = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/DNA_annotation.tsv',
                           header = F, sep = '\t')

hist_expansion = ggplot(genesOfInterest_DNA, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/DNA_repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest_DNA <- as.character(unique(genesOfInterest_DNA$V6))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList_DNA <- factor(as.integer(geneUniverse %in% genesOfInterest_DNA))
names(geneList_DNA) <- geneUniverse

myGOdata_DNA <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_DNA,  annot = annFUN.gene2GO, gene2GO = geneID2GO)

test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_DNA <- getSigGroups(myGOdata_DNA, test.stat)

pvalFis_DNA = score(resultFisher_DNA)
hist(pvalFis_DNA, xlab = "p-values")

allRes_DNA <- GenTable(myGOdata_DNA, classicFisher = resultFisher_DNA, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
str(allRes_DNA)
allRes_DNA$classicFisher = as.numeric(allRes_DNA$classicFisher)

plot_go_DNA = ggplot(allRes_DNA, aes(reorder(Term, Significant), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO enrichment analysis of genes closest to DNA repeats in Hydra vulgaris genome',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_DNA$Significant)) +
  ylim(0,20)


ggsave(plot = plot_go_DNA, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/DNA_transposons_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')

#------------------------------LINE transposons--------------------------------------------------------------------

genesOfInterest_LINE = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_annotation.tsv',
                           header = F, sep = '\t')

hist_expansion_LINE = ggplot(genesOfInterest_LINE, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion_LINE, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest_LINE <- as.character(unique(genesOfInterest_LINE$V6))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList_LINE <- factor(as.integer(geneUniverse %in% genesOfInterest_LINE))
names(geneList_LINE) <- geneUniverse

myGOdata_LINE <- new("topGOdata", description="My project", ontology="MF",
                     allGenes=geneList_LINE,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


resultFisher_LINE <- getSigGroups(myGOdata_LINE, test.stat)
pvalFis_LINE = score(resultFisher_LINE)
hist(pvalFis_LINE, xlab = "p-values")

allRes_LINE <- GenTable(myGOdata_LINE, classicFisher = resultFisher_LINE,
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
str(allRes_LINE)
allRes_LINE$classicFisher = as.numeric(allRes_LINE$classicFisher)

plot_go_LINE = ggplot(allRes_LINE, aes(reorder(Term, Significant), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO enrichment analysis of genes closest to LINE repeats in Hydra vulgaris genome',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_LINE$Significant)) +
  ylim(0,20)


ggsave(plot = plot_go_LINE, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/LINE_transposons_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')


#--------------LTR repeats-------------------------------------------------------

genesOfInterest_LTR = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LTR_annotation.tsv',
                                header = F, sep = '\t')

hist_expansion_LTR = ggplot(genesOfInterest_LTR, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion_LTR, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LTR_repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest_LTR <- as.character(unique(genesOfInterest_LTR$V6))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList_LTR <- factor(as.integer(geneUniverse %in% genesOfInterest_LTR))
names(geneList_LTR) <- geneUniverse

myGOdata_LTR <- new("topGOdata", description="My project", ontology="MF",
                     allGenes=geneList_LTR,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


resultFisher_LTR <- getSigGroups(myGOdata_LTR, test.stat)
pvalFis_LTR = score(resultFisher_LTR)
hist(pvalFis_LTR, xlab = "p-values")

allRes_LTR <- GenTable(myGOdata_LTR, classicFisher = resultFisher_LTR,
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
str(allRes_LTR)
allRes_LTR$classicFisher = as.numeric(allRes_LTR$classicFisher)

plot_go_LTR = ggplot(allRes_LTR, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO enrichment analysis of genes closest to LTR repeats in Hydra vulgaris genome',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_LTR$Significant)) +
  ylim(0,20)


ggsave(plot = plot_go_LTR, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/LTR_transposons_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')

#-----------------------RC repeats---------------------------------------------

genesOfInterest_RC = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/RC_annotation.tsv',
                               header = F, sep = '\t')

hist_expansion_RC = ggplot(genesOfInterest_RC, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion_RC, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/RC_repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest_RC <- as.character(unique(genesOfInterest_RC$V6))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList_RC <- factor(as.integer(geneUniverse %in% genesOfInterest_RC))
names(geneList_RC) <- geneUniverse

myGOdata_RC <- new("topGOdata", description="My project", ontology="MF",
                    allGenes=geneList_RC,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


resultFisher_RC <- getSigGroups(myGOdata_RC, test.stat)
pvalFis_RC = score(resultFisher_RC)
hist(pvalFis_RC, xlab = "p-values")

allRes_RC <- GenTable(myGOdata_RC, classicFisher = resultFisher_RC,
                       orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
str(allRes_RC)
allRes_RC$classicFisher = as.numeric(allRes_RC$classicFisher)

plot_go_RC = ggplot(allRes_RC, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO enrichment analysis of genes closest to RC repeats in Hydra vulgaris genome',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_RC$Significant)) +
  ylim(0,20)


ggsave(plot = plot_go_RC, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/RC_transposons_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')


 #-----------------------Other transposons plot------------------

genesOfInterest_Other = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/Other_annotation.tsv',
                              header = F, sep = '\t')

hist_expansion_Other = ggplot(genesOfInterest_Other, aes(V1)) +
  geom_bar(stat = 'count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  labs(x = 'Contigs', y = 'Number of repeats')

ggsave(plot = hist_expansion_Other, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/Other_repeat_hist.tiff',
       width = 10, height = 10, device = 'tiff')



genesOfInterest_Other <- as.character(unique(genesOfInterest_Other$V6))

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

geneList_Other <- factor(as.integer(geneUniverse %in% genesOfInterest_Other))
names(geneList_Other) <- geneUniverse

myGOdata_Other <- new("topGOdata", description="My project", ontology="MF",
                   allGenes=geneList_Other,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


resultFisher_Other <- getSigGroups(myGOdata_Other, test.stat)
pvalFis_Other = score(resultFisher_Other)
hist(pvalFis_Other, xlab = "p-values")

allRes_Other <- GenTable(myGOdata_Other, classicFisher = resultFisher_Other,
                      orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
str(allRes_Other)
allRes_Other$classicFisher = as.numeric(allRes_Other$classicFisher)

plot_go_Other = ggplot(allRes_Other, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = 'bold'),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size=14,color = "black"),
        axis.title.y = element_text(size=14,color = "black"),
        strip.text.x = element_text(size = 14,colour = "black"),
        legend.text = element_text(size=14,color = "black"),
        legend.title = element_text(size = 14,color = "black")) +
  labs(title = 'GO enrichment analysis of genes closest to Other repeats in Hydra vulgaris genome',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_Other$Significant)) +
  ylim(0,20)


ggsave(plot = plot_go_Other, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/Other_transposons_plot_go.tiff',
       width = 10, height = 10, device = 'tiff')
