#library(dplyr)
library(ggplot2)
library(topGO)
library(RColorBrewer)


#---------readind and reformating dataframe-------
green_hydra = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_repeats_near_genes_filtered.tsv',
                             sep = '\t', stringsAsFactors = F, header = F)

colnames(green_hydra)[20] = 'Gene'
colnames(green_hydra)[19] = 'Distance'
cr1 = green_hydra[grepl('LINE/CR1', green_hydra$V3), ]
cr1$Distance = as.numeric(as.character(cr1$Distance))

#-------------------------Gene ontology annotaton file----------------------------------------------------------
geneID2GO = readMappings(file = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/functional_annotation/GO_annotation_interpro.tsv')

geneUniverse <- names(geneID2GO)
geneUniverse = gsub('\"', '', geneUniverse)
geneUniverse = gsub('"', '', geneUniverse)

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis = function(genesOfInterest) {
  genesOfInterest <- as.character(unique(genesOfInterest$Gene))
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(myGOdata, test.stat)
  pvalFis = score(resultFisher)
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
  allRes$classicFisher = as.numeric(allRes$classicFisher)
  return(allRes)
}

res = GO_analysis(repeats_and_genes)

go_all = ggplot(res, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher),
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
  labs(title = 'Gene ontology of Hydra viridissima genes nearest to repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(res$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 6))

ggsave(plot = plot_go_up_expan_3, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/plot_go_up_expan_3.tiff',
       width = 10, height = 10, device = 'tiff')  

# CR1 repeats

res = GO_analysis(cr1)

res$Term[res$Term == 'apolipoprotein B mRNA editing enzyme com...'] = 'apolipoprotein B mRNA editing enzyme complex'
res$Term[res$Term == 'plasma membrane bounded cell projection ...'] = 'plasma membrane bounded cell projection cytoplasm'


go_all = ggplot(res, aes(reorder(Term, Significant), -log10(classicFisher),
                         fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray8', fill = brewer.pal(11, name = "RdYlGn")[10]) +
  geom_hline(yintercept = c(-log10(0.05)),linetype="dashed") +
  geom_hline(yintercept = c(-log10(0.01)),linetype="dashed") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black"),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size = 18,color = "black")) +
  labs(title = 'Gene ontology of Hydra viridissima genes nearest to CR1 repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 6, 1), limits=c(0, 7))

ggsave(plot = go_all, filename = '/Users/alenakizenko/Documents/PhD/pictures/cr1_go_cc_hvir.tiff',
       width = 10, height = 10, device = 'tiff')  
