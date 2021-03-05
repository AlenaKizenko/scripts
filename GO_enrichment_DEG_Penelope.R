library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(pheatmap)


#-------------------------Gene ontology annotaton file----------------------------------------------------------
geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/GO_annotation/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis = function(df_path) {
  genesOfInterest = read.csv(df_path,
                             header = F, sep = '\t')
  genesOfInterest <- as.character(unique(genesOfInterest$V6))
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(myGOdata, test.stat)
  pvalFis = score(resultFisher)
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
  allRes$classicFisher = as.numeric(allRes$classicFisher)
  return(allRes)
}


#--------------------------DEG loading--------------------------------------------------------------------------

up = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/de_up_gene_ids.txt',
              stringsAsFactors = F, header = F)
colnames(up) = 'Gene'
down = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/de_down_gene_ids.txt',
                stringsAsFactors = F, header = F)

colnames(down) = 'Gene'

#-------------------------Penelope gene annotation loading------------------------------------------------------------

penelope = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/Penelope_gene_annotation.tsv',
                    stringsAsFactors = F, header = F, sep = '\t')
names(penelope)[names(penelope) == 'V19'] <- 'Gene'

#-----------------intersection of UP DEG and Penelope--------------------------------------------------------


penelope_up <- merge(up, penelope, by='Gene')

#-----------------intersection of DOWN DEG and Penelope--------------------------------------------------------


penelope_down <- merge(down, penelope, by='Gene')

#-------normalizing amount of Penelope repeats nearby UP and DOWN DEG----------------------------------------------

de_norm_rep = data.frame('DE_type' = c('Up', 'Down'),
                         'Repeats_num_near_deg' = c(length(penelope_up$Gene), length(penelope_down$Gene)))

deg = ggplot(de_norm_rep, aes(DE_type, Repeats_num_near_deg, fill = DE_type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none') +
  labs(title = 'Amount of Penelope repeats nearby DEG',
       x = '', y = 'Number of repeats') +
  scale_fill_manual(values = brewer.pal(3, 'Paired')[2:1])

ggsave(plot = deg, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/deg_Penelope.tiff',
       width = 10, height = 10, device = 'tiff')

#-------------------UP DEG + Penelope---------------------------

genesOfInterest_up <- as.character(unique(penelope_up$Gene))
geneList_up <- factor(as.integer(geneUniverse %in% genesOfInterest_up))
names(geneList_up) <- geneUniverse
myGOdata_up <- new("topGOdata", description="My project", ontology="MF",
                         allGenes=geneList_up,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_up <- getSigGroups(myGOdata_up, test.stat)
pvalFis_up = score(resultFisher_up)
allRes_up <- GenTable(myGOdata_up, classicFisher = resultFisher_up,
                            orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_up$classicFisher = as.numeric(allRes_up$classicFisher)

allRes_up$Term[allRes_up$Term == 'L-tryptophan transmembrane transporter a...'] = 'L-tryptophan transmembrane transporter activity'
allRes_up$Term[allRes_up$Term == 'N-acylsphingosine amidohydrolase activit...'] = 'N-acylsphingosine amidohydrolase activity'
allRes_up$Term[allRes_up$Term == 'alanine transmembrane transporter activi...'] = 'alanine transmembrane transporter activity'
allRes_up$Term[allRes_up$Term == 'L-proline transmembrane transporter acti...'] = 'L-proline transmembrane transporter activity'
allRes_up$Term[allRes_up$Term == 'hydrolase activity, acting on glycosyl b...'] = 'hydrolase activity, acting on glycosyl bonds'



plot_go_up = ggplot(allRes_up, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'GOs of UP DEG nearby Penelope repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_up$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,5)

ggsave(plot = plot_go_up, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/deg_Penelope_up.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------------DOWN DEG + Penelope---------------------------

genesOfInterest_down <- as.character(unique(penelope_down$Gene))
geneList_down <- factor(as.integer(geneUniverse %in% genesOfInterest_down))
names(geneList_down) <- geneUniverse
myGOdata_down <- new("topGOdata", description="My project", ontology="MF",
                   allGenes=geneList_down,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_down <- getSigGroups(myGOdata_down, test.stat)
pvalFis_down = score(resultFisher_down)
allRes_down <- GenTable(myGOdata_down, classicFisher = resultFisher_down,
                      orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_down$classicFisher = as.numeric(allRes_down$classicFisher)

allRes_down$Term[allRes_down$Term == 'purine phosphoribosyltransferase activit...'] = 'purine phosphoribosyltransferase activity'
allRes_down$Term[allRes_down$Term == 'cyclic-nucleotide phosphodiesterase acti...'] = 'cyclic-nucleotide phosphodiesterase activity'
allRes_down$Term[allRes_down$Term == 'transferase activity, transferring phosp...'] = 'phosphorus groups transferase activity'


plot_go_down = ggplot(allRes_down, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'GOs of DOWN DEG nearby Penelope repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_down$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,5)

ggsave(plot = plot_go_down, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/deg_Penelope_down.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------heatmap of DEG and Penelope repeats----------------

a = full_join(allRes_up[, c(2,6)], allRes_down[, c(2,6)],
              by = 'Term', suffix = c('.up', '.down'))

rownames(a) = a$Term
a = a[, -1]

a_log10 = as.data.frame(sapply(a, FUN = function(x) {-log10(x)}))
rownames(a_log10) = rownames(a)
colnames(a_log10) = c("UP", "DOWN")

a_log10[is.na(a_log10)] = 0
penelope_deg = pheatmap(a_log10, cluster_rows = F, cluster_cols = T,
                   color = c('gray88', brewer.pal(9, "YlOrBr")[2:9]),
                   angle_col = 45, fontsize = 11, fontsize_row = 12, fontsize_col = 12,
                   cellwidth = 60, cellheight = 20,
                   legend_breaks = c(0, 2, 4, max(a_log10)),
                   legend_labels = c("0", "2", "4", "-log10(p-values)\n"))
dev.off()

ggsave(plot = penelope_deg, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/Penelope_heatmap_mf_DEG.tiff',
       width = 10, height = 10, device = 'tiff')

#-------------heatmap of DEG and Penelope repeats: protein number-------------------

a = full_join(allRes_up[, c(2,4)], allRes_down[, c(2,4)],
              by = 'Term', suffix = c('.up', '.down'))

rownames(a) = a$Term
a = a[, -1]

a_log2 = as.data.frame(sapply(a, FUN = function(x) {log2(x)}))
rownames(a_log2) = rownames(a)
colnames(a_log2) = c("UP", "DOWN")

a_log2[is.na(a_log2)] = 0
penelope_deg = pheatmap(a_log2, cluster_rows = F, cluster_cols = T,
                        color = c('gray88', brewer.pal(9, "YlOrBr")[2:9]),
                        angle_col = 45, fontsize = 11, fontsize_row = 12, fontsize_col = 12,
                        cellwidth = 60, cellheight = 20,
                        legend_breaks = c(0, 2, 4, 6, 8, max(a_log2)),
                        legend_labels = c("0", "2", "4", "6", "8", "log2(number of proteins)\n")
                      )
dev.off()

ggsave(plot = penelope_deg, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_Penelope/Penelope_heatmap_mf_DEG_protein.tiff',
       width = 10, height = 10, device = 'tiff')
