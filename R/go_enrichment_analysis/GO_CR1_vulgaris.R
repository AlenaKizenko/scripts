library(ggplot2)
library(topGO)
library(RColorBrewer)


#------CR1 Hydra vulgaris-----------------------

geneID2GO = readMappings(file = '/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/GO_annotation/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

#------------------------Gene Ontology mapping------------------------------------------------------------------

#GO_analysis = function(genesOfInterest) {
#  genesOfInterest <- as.character(unique(genesOfInterest$Gene))
 # geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
#  names(geneList) <- geneUniverse
#  myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
 # test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
#  resultFisher <- getSigGroups(myGOdata, test.stat)
 # pvalFis = score(resultFisher)
  #allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
#  allRes$classicFisher = as.numeric(allRes$classicFisher)
 # return(allRes)
#}

df_brown = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/genome_annotation_all_repeats/first_round/LINE_CR1_fr_genes_closest_filt.tab',
                    sep = '\t', header = F)

cr1 = df_brown[grepl('LINE/CR1', df_brown$V3), ]
hist(cr1$V19)

#-----------cr1 all genes analysis---------------
colnames(cr1)[20] = 'Gene'
#res_cr1 = GO_analysis(cr1)

genesOfInterest <- as.character(unique(cr1$Gene))
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(myGOdata, test.stat)
pvalFis = score(resultFisher)
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes$classicFisher = as.numeric(allRes$classicFisher)

allRes$Term[allRes$Term == 'regulation of cellular component movemen...'] = 'regulation of cellular component movement'
allRes$Term[allRes$Term == 'cellular component organization or bioge...'] = 'cellular component organization or biogenesis'
allRes$Term[allRes$Term == 'modification-dependent macromolecule cat...'] = 'modification-dependent macromolecule catabolic process'


go_cr1_hvul = ggplot(allRes, aes(reorder(Term, Significant), -log10(classicFisher),
                                                 fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8', alpha = 0.8, fill = brewer.pal(11, 'Spectral')[10]) +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,color = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = '',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1)
 # scale_y_continuous(breaks=seq(0, 6, 1), limits=c(0, 7))
go_cr1_hvul
#ggsave(plot = go_cr1_hvul, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/go_cr1_bp_hvul.tiff',
 #      device = 'tiff', width = 10, height = 10)

ann.genes <- genesInTerm(myGOdata)

allRes$GO.ID

#fisher.go <- names(score(resultFisher))[1:5]
movement1 <- genesInTerm(myGOdata, whichGO='GO:0051270')[[1]]
movement2 <- genesInTerm(myGOdata, whichGO='GO:0040012')[[1]]
mov = intersect(movement1, movement2)

adhesion1 <- genesInTerm(myGOdata, whichGO='GO:0007155')[[1]]
adhesion2 <- genesInTerm(myGOdata, whichGO='GO:0022610')[[1]]
adh = intersect(adhesion1, adhesion2)

cr1_te_count = as.data.frame(table(cr1$Gene))
cr1_te_count_sub = subset(cr1_te_count, Freq >= 50)

intersect(cr1_te_count_sub$Var1, adh)
intersect(cr1_te_count_sub$Var1, mov)

test = subset(cr1, Gene == 'Sc4wPfr_108.g783.t1')

