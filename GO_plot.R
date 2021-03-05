library(dplyr)
library(ggplot2)
library(topGO)
library(RColorBrewer)


#-------------------------Gene ontology annotaton file----------------------------------------------------------
geneID2GO = readMappings(file = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/GO_annotation/reformatted_go.csv')

geneUniverse <- names(geneID2GO) 

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis = function(genesOfInterest) {
  genesOfInterest <- as.character(unique(genesOfInterest$Gene))
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


#------3 hours stage GO analysis-----------

# GO 3h expan up
GO_3_up_expan = GO_analysis(expan_up_3)

plot_go_up_expan_3 = ggplot(GO_3_up_expan, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher),
                                               fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8') +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.75, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Gene ontology of upregulated (3h) genes nearby expanded repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(GO_3_up_expan$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4))

ggsave(plot = plot_go_up_expan_3, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/plot_go_up_expan_3.tiff',
       width = 10, height = 10, device = 'tiff')  