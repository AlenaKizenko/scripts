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

DNA = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/DNA_annotation.tsv')
LINE = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/LINE_annotation.tsv')
LTR = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/LTR_annotation.tsv')
Other = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/Other_annotation.tsv')
RC = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/RC_annotation.tsv')
Unknown = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/repeat_class_annotation/Unknown_annotation.tsv')



a = full_join(DNA[, c(2,6)], LINE[, c(2,6)],
              by = 'Term', suffix = c('.DNA', '.LINE'))
b = full_join(LTR[, c(2,6)], Other[, c(2,6)],
              by = 'Term', suffix = c('.LTR', '.Other'))
c = full_join(RC[, c(2,6)], Unknown[, c(2,6)],
              by = 'Term', suffix = c('.RC', '.Unknown'))

ab = full_join(a, b, by = 'Term', suffix = c('', ''))
abc = full_join(ab, c, by = 'Term', suffix = c('', ''))


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


rownames(abc) = abc$Term
abc = abc[, -1]

abc_log10 = as.data.frame(sapply(abc, FUN = function(x) {-log10(x)}))
rownames(abc_log10) = rownames(abc)
colnames(abc_log10) = c("DNA", "LINE", "LTR",
                        "Other", "RC", "Unknown")
#abcde_log10 = as.data.frame(t(abcde_log10))


abc_log10[is.na(abc_log10)] = 0
lines_types = pheatmap(abc_log10, cluster_rows = T, cluster_cols = T,
                       color = colorRampPalette(brewer.pal(9, "Purples"))(50),
                       na_col = 'gray77', angle_col = 45, fontsize_row = 12, fontsize_col = 10)

ggsave(plot = lines_types, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/LINE_types_annotation/lines_types.tiff',
       width = 10, height = 10, device = 'tiff')

#----------without other and unknown---------------------------------------


d = full_join(DNA[, c(2,6)], LINE[, c(2,6)],
              by = 'Term', suffix = c('.DNA', '.LINE'))
e = full_join(LTR[, c(2,6)], RC[, c(2,6)],
              by = 'Term', suffix = c('.LTR', '.RC'))

de = full_join(d, e, by = 'Term', suffix = c('', ''))

de$Term[de$Term == 'transferase activity, transferring acyl ...'] = 'transferase activity, transferring acyl groups'
de$Term[de$Term == 'transferase activity, transferring phosp...'] = 'transferase activity, transferring phosphorus-containing groups'
de$Term[de$Term == 'transcription regulatory region DNA bind...'] = 'transcription regulatory region DNA binding'
de$Term[de$Term == 'transmembrane receptor protein phosphata...'] = 'transmembrane receptor protein phosphatase activity'


rownames(de) = de$Term
de = de[, -1]
str(de)
de_log10 = as.data.frame(sapply(de, FUN = function(x) {-log10(x)}))
rownames(de_log10) = rownames(de)
colnames(de_log10) = c("DNA", "LINE", "LTR", "RC")
#abcde_log10 = as.data.frame(t(abcde_log10))


de_log10[is.na(de_log10)] = 0
lines_types = pheatmap(de_log10, cluster_rows = T, cluster_cols = T, cellwidth = 20, cellheight = 20,
                       color = c('gray88', brewer.pal(9, "Purples")[2:9]),
                       na_col = 'gray77', angle_col = 45, fontsize_row = 12, fontsize_col = 12)

ggsave(plot = lines_types, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_all_repeats/classes_annotation.tiff',
       width = 10, height = 10, device = 'tiff')
