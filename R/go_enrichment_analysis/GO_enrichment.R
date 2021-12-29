library(topGO)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(ggsignif)


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

#-------------------------expansion------------------------------------------------------------

exp = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/expansion/repeat_expansion_annotation_selected.gff')
str(exp)

exp$Term[exp$Term == 'oxidoreductase activity, acting on CH-OH...'] = 'oxidoreductase activity'
exp$Term[exp$Term == 'cellular component organization or bioge...'] = 'cellular component organization or biogenesis'
exp$Term[exp$Term == 'regulation of cellular component movemen...'] = 'regulation of cellular component movement'
exp$Term[exp$Term == 'movement of cell or subcellular componen...'] = 'movement of cell or subcellular component'


plot_go = ggplot(exp, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Repeat expansion genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(exp$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,6)


ggsave(plot = plot_go, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/expansion_plot_go_cc.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------------------rest------------------------------------------------------------


rest = GO_analysis('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/common/rest_repeats_annotation_selected.gff')


rest$Term[rest$Term == 'transferase activity, transferring phosp...'] = 'transferase activity'
rest$Term[rest$Term == 'transmembrane receptor protein kinase ac...'] = 'transmembrane receptor protein kinase activity'
rest$Term[rest$Term == 'calcium ion transmembrane transporter ac...'] = 'calcium ion transmembrane transporter activity'

rest$Term[rest$Term == 'movement of cell or subcellular componen...'] = 'movement of cell or subcellular component'


plot_go_rest = ggplot(rest, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Rest repeat genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(rest$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
 ylim(0,8)


ggsave(plot = plot_go_rest, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/rest_plot_go_cc.tiff',
       width = 10, height = 10, device = 'tiff')

#--------------merged GO------------------------------

a = full_join(exp[, c(2,6)], rest[, c(2,6)],
              by = 'Term', suffix = c('.exp', '.rest'))

a$Term[a$Term == 'oxidoreductase activity, acting on CH-OH...'] = 'oxidoreductase activity, acting on CH-OH group of donors'
a$Term[a$Term == 'transferase activity, transferring phosp...'] = 'transferase activity, transferring phosphorus-containing groups'
a$Term[a$Term == 'transmembrane receptor protein kinase ac...	'] = 'transmembrane receptor protein kinase activity'
a$Term[a$Term == 'calcium ion transmembrane transporter ac...'] = 'calcium ion transmembrane transporter activity'

rownames(a) = a$Term
a = a[, -1]

a_log10 = as.data.frame(sapply(a, FUN = function(x) {-log10(x)}))
rownames(a_log10) = rownames(a)
colnames(a_log10) = c("Expanded", "Non-expanded")

a_log10[is.na(a_log10)] = 0
cr1 = pheatmap(a_log10, cluster_rows = T, cluster_cols = T,
                       color = c('gray88', brewer.pal(9, "Blues")[2:9]),
                       angle_col = 45, fontsize = 12, fontsize_row = 14, fontsize_col = 14,
                      legend_breaks = c(0, 2, 4, 6, 8, 10, 12, 14, max(a_log10)),
                      legend_labels = c("0", "2", "4", "6", "8", "10", "12", "14", "-log10(p-values)\n"))
#dev.off()

ggsave(plot = cr1, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/CR1_heatmap_bp.tiff',
       width = 10, height = 10, device = 'tiff')


  #------------------test for wnt genes-------------------------------------------------------------------

wnt = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/wnt_genes.csv', stringsAsFactors = F)



expan = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_repeats_annotation_expansion.tsv',
                 stringsAsFactors = F, header = F, sep = '\t')
rst = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_repeats_annotation_rest.tsv',
               stringsAsFactors = F, header = F, sep = '\t')

names(expan)[names(expan) == 'V19'] <- 'Gene'
names(rst)[names(rst) == 'V19'] <- 'Gene'

#colnames(expan) = c('Chromosome_rep', 'Repeat', 'Start_rep',
 #                   'End_rep', 'Chromosome_gene', 'Gene', 'Start_gene',
  #                  'End_gene', 'Distance')

#colnames(rst) = c('Chromosome_rep', 'Repeat', 'Start_rep',
 #       'End_rep', 'Chromosome_gene', 'Gene', 'Start_gene',
  #      'End_gene', 'Distance')

colnames(wnt) = c("NCBI_ID", "Gene_name", "Gene")

expan_wnt <- merge(wnt, expan, by='Gene')

rst_wnt <- merge(wnt, rst, by='Gene')

  expan_genes = as.data.frame(unique(expan_wnt$Gene))
  wnt_genes = as.data.frame(unique(wnt$Gene))
  rst_genes = as.data.frame(unique(rst_wnt$Gene))

df = data.frame('Expansion' = c(6, 7), 'Rest' = c(0, 13))
rownames(df) = c('Near repeats', 'Away from repeats')
fisher.test(df)

# DE genes analysis

expan = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_repeats_annotation_expansion.tsv',
                 stringsAsFactors = F, header = F, sep = '\t')
rst = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_repeats_annotation_rest.tsv',
               stringsAsFactors = F, header = F, sep = '\t')

names(expan)[names(expan) == 'V19'] <- 'Gene'
names(rst)[names(rst) == 'V19'] <- 'Gene'

up = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/de_up_gene_ids.txt',
              stringsAsFactors = F, header = F)
colnames(up) = 'Gene'
down = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/de_down_gene_ids.txt',
                stringsAsFactors = F, header = F)
colnames(down) = 'Gene'
random = read.csv('/home/alena/Documents/PhD/hydra_vulgaris_chrom/genes_headers_4perc.txt',
                  header = F, stringsAsFactors = F)
random <- as.data.frame(unlist(lapply(random, gsub, pattern=">", replacement="")))
colnames(random) = 'Gene'

expan_up <- merge(up, expan, by='Gene')
expan_down <- merge(down, expan, by='Gene')
expan_random = merge(random, expan, by='Gene')

#expan_down_unique = unique(expan_down$Gene)
#write.csv(expan_down_unique, '/home/alena/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/expan_down_genes.txt')

rst_up <- merge(up, rst, by='Gene')
rst_down <- merge(down, rst, by='Gene')
rst_random <- merge(random, rst, by='Gene')

rep_genes_fisher = data.frame('Up' = c(length(expan_up$Gene), length(rst_up$Gene)),
                              'Down' = c(length(expan_down$Gene), length(rst_down$Gene)),
                              'Random' = c(length(expan_random$Gene), length(rst_random$Gene)))
rownames(rep_genes_fisher) = c('Expansion', 'Rest')

fish = fisher.test(rep_genes_fisher[,c(2,3)], alternative = 'two.sided')
fish$p.value

View(head(expan))

de_norm_rep = data.frame('DE_type' = c('Upregulated genes', 'Downregulated genes', 'Random genes', 'Upregulated genes', 'Downregulated genes', 'Random genes'),
  'Repeats_num_near_deg' = c(length(expan_up$Gene)/length(expan$V1), length(expan_down$Gene)/length(expan$V1), length(expan_random$Gene)/length(expan$V1),
           length(rst_up$Gene)/length(rst$V1), length(rst_down$Gene)/length(rst$V1), length(rst_random$Gene)/length(rst$V1)),
  'Repeat_type' = c('Expansion', 'Expansion', 'Expansion', 'Rest', 'Rest', 'Rest'),
  'Order' = c(3,1,5,4,2,6))

deg = ggplot(de_norm_rep, aes(reorder(DE_type, Order), Repeats_num_near_deg*100, fill = Repeat_type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=16,color = "black"),
        axis.text.y = element_text(size=16,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'CR1 repeats nearby differentially expressed genes during regeneration',
       x = '', y = 'Percentage of repeats',
       fill = 'Repeat type') +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6")) +
  geom_text(x=1, y=6.8, label=paste("p-value=",round(fisher.test(rep_genes_fisher[,c(2,3)],
                                                                 alternative = 'greater')$p.value, 6)), size = 6, color = 'gray7') +
  geom_text(x=2, y=3.8, label=paste("p-value=",round(fisher.test(rep_genes_fisher[,c(1,3)],
                                                                 alternative = 'two.sided')$p.value, 6)), size = 6, color = 'gray7') +
  
  geom_text(x=3, y=6.8, label="One-sided Fisher exact test", size = 5, color = 'gray7')    


ggsave(plot = deg, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/deg_and_random_CR1.tiff',
       width = 10, height = 10, device = 'tiff')

# DEG analysis GO

geneID2GO = readMappings(file = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

#-------------expansion UP DEG analysis--------------------------------

genesOfInterest_expan_up <- as.character(unique(expan_up$Gene))
geneList_expan_up <- factor(as.integer(geneUniverse %in% genesOfInterest_expan_up))
names(geneList_expan_up) <- geneUniverse
myGOdata_expan_up <- new("topGOdata", description="My project", ontology="MF",
                         allGenes=geneList_expan_up,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_expan_up <- getSigGroups(myGOdata_expan_up, test.stat)
pvalFis_expan_up = score(resultFisher_expan_up)
allRes_expan_up <- GenTable(myGOdata_expan_up, classicFisher = resultFisher_expan_up,
                            orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_expan_up$classicFisher = as.numeric(allRes_expan_up$classicFisher)

allRes_expan_up$Term[allRes_expan_up$Term == 'transferase activity, transferring sulfu...'] = 'transferase activity, transferring sulfur-containing groups'
allRes_expan_up$Term[allRes_expan_up$Term == 'transferase activity, transferring glyco...'] = 'transferase activity, transferring glycosyl groups'
allRes_expan_up$Term[allRes_expan_up$Term == 'carboxylic acid transmembrane transporte...'] = 'carboxylic acid transmembrane transporter activity'
allRes_expan_up$Term[allRes_expan_up$Term == 'N-acylsphingosine amidohydrolase activit...'] = 'N-acylsphingosine amidohydrolase activity'
allRes_expan_up$Term[allRes_expan_up$Term == 'alanine transmembrane transporter activi...'] = 'alanine transmembrane transporter activity'



plot_go_up_expan = ggplot(allRes_expan_up, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Repeat expansion UP genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_expan_up$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,10)

ggsave(plot = plot_go_up_expan, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/deg_expansion_up.tiff',
       width = 10, height = 10, device = 'tiff')

#--------------rest UP DEG analysis---------------------

genesOfInterest_rest_up <- as.character(unique(rst_up$Gene))
geneList_rest_up <- factor(as.integer(geneUniverse %in% genesOfInterest_rest_up))
names(geneList_rest_up) <- geneUniverse
myGOdata_rest_up <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_rest_up,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_rest_up <- getSigGroups(myGOdata_rest_up, test.stat)
pvalFis_rest_up = score(resultFisher_rest_up)
allRes_rest_up <- GenTable(myGOdata_rest_up, classicFisher = resultFisher_rest_up,
                           orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_rest_up$classicFisher = as.numeric(allRes_rest_up$classicFisher)


allRes_rest_up$Term[allRes_rest_up$Term == 'beta-ketoacyl-acyl-carrier-protein synth...'] = 'beta-ketoacyl-acyl-carrier-protein synthase II activity'
allRes_rest_up$Term[allRes_rest_up$Term == 'ammonium transmembrane transporter activ...'] = 'ammonium transmembrane transporter activity'
allRes_rest_up$Term[allRes_rest_up$Term == 'carbon dioxide transmembrane transporter...'] = 'carbon dioxide transmembrane transporter activity'


plot_go_up_rest = ggplot(allRes_rest_up, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Rest repeat UP genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_rest_up$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,10)


ggsave(plot = plot_go_up_rest, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/deg_rest_up.tiff',
       width = 10, height = 10, device = 'tiff')

#-------expansion DOWN DEG analysis-------------------------------

genesOfInterest_expan_down <- as.character(unique(expan_down$Gene))
geneList_expan_down <- factor(as.integer(geneUniverse %in% genesOfInterest_expan_down))
names(geneList_expan_down) <- geneUniverse
myGOdata_expan_down <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_expan_down,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_expan_down <- getSigGroups(myGOdata_expan_down, test.stat)
pvalFis_expan_down = score(resultFisher_expan_down)
allRes_expan_down <- GenTable(myGOdata_expan_down, classicFisher = resultFisher_expan_down,
                              orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_expan_down$classicFisher = as.numeric(allRes_expan_down$classicFisher)

allRes_expan_down$Term[allRes_expan_down$Term == 'cyclic nucleotide-dependent protein kina...'] = 'cyclic nucleotide-dependent protein kinase activity'


plot_go_down_expan = ggplot(allRes_expan_down, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Repeat expansion DOWN genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_expan_down$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,10)

ggsave(plot = plot_go_down_expan, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/deg_expansion_down.tiff',
       width = 10, height = 10, device = 'tiff')

#-----------rest DOWN DEG analysis-------------------------------

genesOfInterest_rest_down <- as.character(unique(rst_down$Gene))
geneList_rest_down <- factor(as.integer(geneUniverse %in% genesOfInterest_rest_down))
names(geneList_rest_down) <- geneUniverse
myGOdata_rest_down <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_rest_down,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_rest_down <- getSigGroups(myGOdata_rest_down, test.stat)
pvalFis_rest_down = score(resultFisher_rest_down)
allRes_rest_down <- GenTable(myGOdata_rest_down, classicFisher = resultFisher_rest_down,
                             orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_rest_down$classicFisher = as.numeric(allRes_rest_down$classicFisher)


allRes_rest_down$Term[allRes_rest_down$Term == 'oxidoreductase activity, acting on the C...'] = 'oxidoreductase activity'
allRes_rest_down$Term[allRes_rest_down$Term == 'alpha-1,6-mannosylglycoprotein 6-beta-N-...'] = '6-beta-N-acetylglucosaminyltransferase activity'



plot_go_down_rest = ggplot(allRes_rest_down, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher), fill = Significant)) +
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
  labs(title = 'Rest repeat DOWN genes\' GOs',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(allRes_rest_down$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  ylim(0,10)


ggsave(plot = plot_go_down_rest, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/deg_rest_down.tiff',
       width = 10, height = 10, device = 'tiff')

#-------------heatmap of DEG and CR1 repeats----------------

a = full_join(allRes_expan_up[, c(2,6)], allRes_expan_down[, c(2,6)],
              by = 'Term', suffix = c('.expup', '.expdown'))

b = full_join(allRes_rest_up[, c(2,6)], allRes_rest_down[, c(2,6)],
              by = 'Term', suffix = c('.restup', '.restdown'))

c = full_join(a, b, by = 'Term', suffix = c('', ''))

rownames(c) = c$Term
c = c[, -1]

c_log10 = as.data.frame(sapply(c, FUN = function(x) {-log10(x)}))
rownames(c_log10) = rownames(c)
colnames(c_log10) = c("Expanded UP", "Expanded DOWN", "Rest UP", "Rest DOWN")

c_log10[is.na(c_log10)] = 0
cr1_deg = pheatmap(c_log10, cluster_rows = F, cluster_cols = F,
               color = c('gray88', brewer.pal(9, "PuRd")[2:9]),
               cellwidth = 40, cellheight = 15,
               angle_col = 45, fontsize = 11, fontsize_row = 12, fontsize_col = 12,
               legend_breaks = c(0, 2, 4, 6, 8, 10, max(c_log10)),
               legend_labels = c("0", "2", "4", "6", "8", "10", "-log10(p-values)\n"))
dev.off()

ggsave(plot = cr1_deg, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/CR1_heatmap_mf_DEG.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------heatmap of DEG and CR1 repeats: protein number-------------------

a = full_join(allRes_expan_up[, c(2,4)], allRes_expan_down[, c(2,4)],
              by = 'Term', suffix = c('.expup', '.expdown'))

b = full_join(allRes_rest_up[, c(2,4)], allRes_rest_down[, c(2,4)],
              by = 'Term', suffix = c('.restup', '.restdown'))

c = full_join(a, b, by = 'Term', suffix = c('', ''))

rownames(c) = c$Term
c = c[, -1]

c_log2 = as.data.frame(sapply(c, FUN = function(x) {log2(x)}))
rownames(c_log2) = rownames(c)
colnames(c_log2) = c("Expanded UP", "Expanded DOWN", "Rest UP", "Rest DOWN")

c_log2[is.na(c_log2)] = 0

cr1_deg_prot = pheatmap(c_log2, cluster_rows = F, cluster_cols = F,
                        color = c('gray88', brewer.pal(9, "PuRd")[2:9]),
                        angle_col = 45, fontsize = 11, fontsize_row = 12, fontsize_col = 12,
                        cellwidth = 40, cellheight = 15,
                        legend_breaks = c(0, 2, 4, 6, 8, max(c_log2)),
                        legend_labels = c("0", "2", "4", "6", "8", "log2(number of proteins)\n")
)
dev.off()

ggsave(plot = cr1_deg_prot, filename = '/home/alena/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/CR1_heatmap_mf_DEG_protein.tiff',
       width = 10, height = 10, device = 'tiff')
