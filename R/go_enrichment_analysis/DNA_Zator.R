library(ggplot2)
library(topGO)
library(RColorBrewer)


#---------
# Green Hydra
green_hydra = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeat_genes_intersection/hvir_repeats_near_genes_filtered.tsv',
                       sep = '\t', stringsAsFactors = F, header = F)

colnames(green_hydra)[20] = 'Gene'
colnames(green_hydra)[19] = 'Distance'

pao_hvir = green_hydra[grepl('DNA/Zator', green_hydra$V3), ]
pao_hvir$Distance = as.numeric(as.character(pao_hvir$Distance))

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis_hvir = function(genesOfInterest) {
  geneID2GO_hvir = readMappings(file = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/functional_annotation/GO_annotation_interpro.tsv')
  geneUniverse_hvir = names(geneID2GO_hvir)
  geneUniverse_hvir = gsub('\"', '', geneUniverse_hvir)
  geneUniverse_hvir = gsub('"', '', geneUniverse_hvir)
  
  genesOfInterest <- as.character(unique(genesOfInterest$Gene))
  geneList <- factor(as.integer(geneUniverse_hvir %in% genesOfInterest))
  names(geneList) <- geneUniverse_hvir
  myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO_hvir)
  test.stat <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(myGOdata, test.stat)
  pvalFis = score(resultFisher)
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
  allRes$classicFisher = as.numeric(allRes$classicFisher)
  return(allRes)
}

res_hvir = GO_analysis_hvir(pao_hvir)

res_hvir$Term[res_hvir$Term == 'cellular aromatic compound metabolic pro...'] = 'cellular aromatic compound metabolic process'
res_hvir$Term[res_hvir$Term == 'thiamine diphosphate biosynthetic proces...'] = 'thiamine diphosphate biosynthetic process'
res_hvir$Term[res_hvir$Term == 'cellular aromatic compound metabolic pro...'] = 'ccellular aromatic compound metabolic process'

go_all_hvir = ggplot(res_hvir, aes(reorder(Term, Significant), -log10(classicFisher),
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
  labs(title = 'Gene ontology of Hydra viridissima genes nearest to DNA/Zator repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(res_hvir$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 10, 1), limits=c(0, 11))

ggsave(plot = go_all_hvir, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/pictures/plot_go_hvir.tiff',
       width = 10, height = 10, device = 'tiff')  

# Brown Hydra

brown_hydra = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1/all_repeats_gene_annotation_filtered_dist.tsv',
                       sep = '\t', stringsAsFactors = F, header = F)

colnames(brown_hydra)[20] = 'Gene'
colnames(brown_hydra)[19] = 'Distance'
pao_brown = brown_hydra[grepl('DNA/Zator', brown_hydra$V3), ]
pao_brown$Distance = as.numeric(as.character(pao_brown$Distance))

#------------------------Gene Ontology mapping------------------------------------------------------------------

GO_analysis = function(genesOfInterest) {
  geneID2GO = readMappings(file = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/GO_annotation/reformatted_go.csv')
  geneUniverse <- names(geneID2GO) 
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

res_hvulg = GO_analysis(pao_brown)

res_hvulg$Term[res_hvulg$Term == 'developmental process involved in reprod...'] = 'developmental process involved in reproduction'

go_all = ggplot(res_hvulg, aes(reorder(Term, Significant), -log10(classicFisher),
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
  labs(title = 'Gene ontology of Hydra vulgaris genes nearest to DNA/Zator repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(res_hvulg$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 10, 1), limits=c(0, 11))

ggsave(plot = go_all, filename = '/Users/alenakizenko/Documents/PhD/CR1_green_brown_comparison/pictures/go_all_hvulg.tiff',
       width = 10, height = 10, device = 'tiff')  
