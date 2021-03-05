library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(pheatmap)


#-------------------------Gene ontology annotaton file----------------------------------------------------------
geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis = function(df_path) {
  genesOfInterest = read.csv(df_path,
                                 header = F, sep = '\t')
  genesOfInterest <- as.character(unique(genesOfInterest$V6))
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(myGOdata, test.stat)
  pvalFis = score(resultFisher)
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
  allRes$classicFisher = as.numeric(allRes$classicFisher)
  return(allRes)
}

CR1 = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_CR1_annotation.tsv')
CR1_Zenon = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_CR1-Zenon_annotation.tsv')
L1 = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_L1_annotation.tsv')
L2 = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_L2_annotation.tsv')
LINE = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_LINE_annotation.tsv')
Penelope = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_Penelope_annotation.tsv')
Proto2 = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_Proto2_annotation.tsv')
R2 = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_R2_annotation.tsv')
Rex_Babar = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_Rex-Babar_annotation.tsv')
RTE_BovB = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/LINE_RTE-BovB_annotation.tsv')


a = full_join(L1[, c(2,6)], L2[, c(2,6)],
              by = 'Term', suffix = c('.L1', '.L2'))
b = full_join(CR1[, c(2,6)], CR1_Zenon[, c(2,6)],
              by = 'Term', suffix = c('.CR1', '.CR1-Zenon'))
c = full_join(LINE[, c(2,6)], Penelope[, c(2,6)],
              by = 'Term', suffix = c('.LINE', '.Penelope'))
d = full_join(Proto2[, c(2,6)], R2[, c(2,6)],
              by = 'Term', suffix = c('.Proto2', '.R2'))
e = full_join(Rex_Babar[, c(2,6)], RTE_BovB[, c(2,6)],
              by = 'Term', suffix = c('.Rex-Babar', '.RTE-BovB'))


ab = full_join(a, b, by = 'Term', suffix = c('', ''))
cd = full_join(c, d, by = 'Term', suffix = c('', ''))
adcd = full_join(ab, cd, by = 'Term', suffix = c('', ''))
abcde = full_join(adcd, e, by = 'Term', suffix = c('', ''))

abcde$Term[abcde$Term == 'ubiquitin-like protein conjugating enzym...'] = 'ubiquitin-like protein conjugating enzyme binding'
abcde$Term[abcde$Term == 'transferase activity, transferring phosp...'] = 'transferase activity, transferring phosphorus-containing groups'
abcde$Term[abcde$Term == 'transcription regulatory region DNA bind...'] = 'transcription regulatory region DNA binding'
abcde$Term[abcde$Term == 'endonuclease activity, active with eithe...'] = 'endonuclease activity, producing 5\'-phosphomonoesters'
#abcde$Term[abcde$Term == 'endonuclease activity, active with eithe...'] = 'endonuclease activity, active with either ribo- or \n deoxyribonucleic acids and producing 5\'-phosphomonoesters'
abcde$Term[abcde$Term == 'protein-containing complex scaffold acti...'] = 'protein-containing complex scaffold activity'
abcde$Term[abcde$Term == 'ubiquitin-like protein transferase activ...'] = 'ubiquitin-like protein transferase activity'
abcde$Term[abcde$Term == 'site-specific endodeoxyribonuclease acti...'] = 'site-specific endodeoxyribonuclease activity, specific for altered base'
abcde$Term[abcde$Term == 'transmembrane receptor protein kinase ac...'] = 'transmembrane receptor protein kinase activity'
abcde$Term[abcde$Term == 'N-acylphosphatidylethanolamine-specific ...'] = 'N-acylphosphatidylethanolamine-specific phospholipase D activity'
abcde$Term[abcde$Term == 'calcium ion transmembrane transporter ac...'] = 'calcium ion transmembrane transporter activity'
abcde$Term[abcde$Term == 'oxidoreductase activity, acting on the C...'] = 'oxidoreductase activity, acting on the CH-NH2 group of donors'
abcde$Term[abcde$Term == 'G-protein beta/gamma-subunit complex bin...'] = 'G-protein beta/gamma-subunit complex binding'


rownames(abcde) = abcde$Term
abcde = abcde[, -1]

abcde_log10 = as.data.frame(sapply(abcde, FUN = function(x) {-log10(x)}))
rownames(abcde_log10) = rownames(abcde)
colnames(abcde_log10) = c("L1", "L2", "CR1", "CR1-Zenon", "LINE", "Penelope",
                          "Proto2", "R2", "Rex-Babar", "RTE-BovB")
#abcde_log10 = as.data.frame(t(abcde_log10))
abcde_log10[is.na(abcde_log10)] = 0

lines_types = pheatmap(abcde_log10, cluster_rows = F, cluster_cols = T,
                     #  cellwidth = 8, cellheight = 8,
                       color = c('gray88', brewer.pal(9, "Greens")[2:9]),
                       na_col = 'gray77', angle_col = 45, fontsize_row = 10, fontsize_col = 10)

ggsave(plot = lines_types, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/lines_types.tiff',
       width = 10, height = 10, device = 'tiff')

