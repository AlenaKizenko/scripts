library(dplyr)
library(ggplot2)
library(topGO)
library(RColorBrewer)

#------------reading parsed data with DEG---------------------
up = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/up_parsed.tsv',
              header = F, stringsAsFactors = F, sep = '\t')

colnames(up) = c('Gene', 'Time', 'DE')

down = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/down_parsed.tsv',
              header = F, stringsAsFactors = F, sep = '\t')
colnames(down) = c('Gene', 'Time', 'DE')

#-------------splitting data by DE-----------------------------
up_true_up = subset(up, DE == 'up')
up_true_down = subset(up, DE == 'down')

down_true_up = subset(down, DE == 'up')
down_true_down = subset(down, DE == 'down')

final_up = rbind(up_true_up, down_true_up)
final_up %>% distinct() -> final_up_unique

final_down = rbind(up_true_down, down_true_down)
final_down %>% distinct() -> final_down_unique

#--------------splitting UP DEG by timecourse-----------------------
str(final_up_unique)

final_up_05 = subset(final_up_unique, Time == 0.5)
final_up_3 = subset(final_up_unique, Time == 3.0)
final_up_6 = subset(final_up_unique, Time == 6.0)
final_up_12 = subset(final_up_unique, Time == 12.0)
final_up_24 = subset(final_up_unique, Time == 24.0)
final_up_48 = subset(final_up_unique, Time == 48.0)

#-------------splitting DOWN DEG by timecourse------------------------

str(final_down_unique)

final_down_05 = subset(final_down_unique, Time == 0.5)
final_down_3 = subset(final_down_unique, Time == 3.0)
final_down_6 = subset(final_down_unique, Time == 6.0)
final_down_12 = subset(final_down_unique, Time == 12.0)
final_down_24 = subset(final_down_unique, Time == 24.0)
final_down_48 = subset(final_down_unique, Time == 48.0)

# LINEs expansion and rest analysis vs UP and DOWN 48h DEG

expan = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_expansion_gene_annotation.tsv',
                 stringsAsFactors = F, header = F, sep = '\t')
View(head(expan))
rst = read.csv('/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/genome_annotation_CR1_new/LINE_CR1_rest_gene_annotation.tsv',
               stringsAsFactors = F, header = F, sep = '\t')

names(expan)[names(expan) == 'V20'] <- 'Gene'
names(rst)[names(rst) == 'V20'] <- 'Gene'

#----------0.5 hours regeneration-----------------
expan_up_05 <- merge(final_up_05, expan, by='Gene')
expan_down_05 <- merge(final_down_05, expan, by='Gene')

rst_up_05 <- merge(final_up_05, rst, by='Gene')
rst_down_05 <- merge(final_down_05, rst, by='Gene')

#----------3 hours regeneration-----------------
expan_up_3 <- merge(final_up_3, expan, by='Gene')
expan_down_3 <- merge(final_down_3, expan, by='Gene')

rst_up_3 <- merge(final_up_3, rst, by='Gene')
rst_down_3 <- merge(final_down_3, rst, by='Gene')

#----------6 hours of regeneration------------
expan_up_6 <- merge(final_up_6, expan, by='Gene')
expan_down_6 <- merge(final_down_6, expan, by='Gene')

rst_up_6 <- merge(final_up_6, rst, by='Gene')
rst_down_6 <- merge(final_down_6, rst, by='Gene')

#----------12 hours regeneration-----------------
expan_up_12 <- merge(final_up_12, expan, by='Gene')
expan_down_12 <- merge(final_down_12, expan, by='Gene')

rst_up_12 <- merge(final_up_12, rst, by='Gene')
rst_down_12 <- merge(final_down_12, rst, by='Gene')

#----------24 hours regeneration-----------------
expan_up_24 <- merge(final_up_24, expan, by='Gene')
expan_down_24 <- merge(final_down_24, expan, by='Gene')

rst_up_24 <- merge(final_up_24, rst, by='Gene')
rst_down_24 <- merge(final_down_24, rst, by='Gene')

#----------48 hours regeneration-----------------
expan_up_48 <- merge(final_up_48, expan, by='Gene')
expan_down_48 <- merge(final_down_48, expan, by='Gene')

rst_up_48 <- merge(final_up_48, rst, by='Gene')
rst_down_48 <- merge(final_down_48, rst, by='Gene')

#----plotting with all stages---------------------

de_norm_rep = data.frame('DE_type' = c('Up 0.5h', 'Down 0.5h', 'Up 0.5h', 'Down 0.5h',
                                       'Up 3h', 'Down 3h', 'Up 3h', 'Down 3h',
                                       'Up 6h', 'Down 6h', 'Up 6h', 'Down 6h',
                                       'Up 12h', 'Down 12h', 'Up 12h', 'Down 12h',
                                       'Up 24h', 'Down 24h', 'Up 24h', 'Down 24h',
                                       'Up 48h', 'Down 48h', 'Up 48h', 'Down 48h'),
                         'Repeats_num_near_deg' = c(length(expan_up_05$Gene)/length(expan$Gene)/length(final_up_05$Gene), length(expan_down_05$Gene)/length(expan$Gene)/length(final_down_05$Gene),
                                                    length(rst_up_05$Gene)/length(rst$Gene)/length(final_up_05$Gene), length(rst_down_05$Gene)/length(rst$Gene)/length(final_down_05$Gene),
                                                    length(expan_up_3$Gene)/length(expan$Gene)/length(final_up_3$Gene), length(expan_down_3$Gene)/length(expan$Gene)/length(final_down_3$Gene),
                                                    length(rst_up_3$Gene)/length(rst$Gene)/length(final_up_3$Gene), length(rst_down_3$Gene)/length(rst$Gene)/length(final_down_3$Gene),
                                                    length(expan_up_6$Gene)/length(expan$Gene)/length(final_up_6$Gene), length(expan_down_6$Gene)/length(expan$Gene)/length(final_down_6$Gene),
                                                    length(rst_up_6$Gene)/length(rst$Gene)/length(final_up_6$Gene), length(rst_down_6$Gene)/length(rst$Gene)/length(final_down_6$Gene),
                                                    length(expan_up_12$Gene)/length(expan$Gene)/length(final_up_12$Gene), length(expan_down_12$Gene)/length(expan$Gene)/length(final_down_12$Gene),
                                                    length(rst_up_12$Gene)/length(rst$Gene)/length(final_up_12$Gene), length(rst_down_12$Gene)/length(rst$Gene)/length(final_down_12$Gene),
                                                    length(expan_up_24$Gene)/length(expan$Gene)/length(final_up_24$Gene), length(expan_down_24$Gene)/length(expan$Gene)/length(final_down_24$Gene),
                                                    length(rst_up_24$Gene)/length(rst$Gene)/length(final_up_24$Gene), length(rst_down_24$Gene)/length(rst$Gene)/length(final_down_24$Gene),
                                                    length(expan_up_48$Gene)/length(expan$Gene)/length(final_up_48$Gene), length(expan_down_48$Gene)/length(expan$Gene)/length(final_down_48$Gene),
                                                    length(rst_up_48$Gene)/length(rst$Gene)/length(final_up_48$Gene), length(rst_down_48$Gene)/length(rst$Gene)/length(final_down_48$Gene)
                                                    ),
                         'Repeat_type' = c('Expansion repeats', 'Expansion repeats', 'Rest repeats',
                                           'Rest repeats', 'Expansion repeats', 'Expansion repeats', 'Rest repeats', 'Rest repeats',
                                           'Expansion repeats', 'Expansion repeats', 'Rest repeats', 'Rest repeats',
                                           'Expansion repeats', 'Expansion repeats', 'Rest repeats', 'Rest repeats',
                                           'Expansion repeats', 'Expansion repeats', 'Rest repeats', 'Rest repeats',
                                           'Expansion repeats', 'Expansion repeats', 'Rest repeats', 'Rest repeats'),
                         'Number_of_repeats' = c(length(expan_up_05$Gene), length(expan_down_05$Gene),
                         length(rst_up_05$Gene), length(rst_down_05$Gene),
                         length(expan_up_3$Gene), length(expan_down_3$Gene),
                         length(rst_up_3$Gene), length(rst_down_3$Gene),
                         length(expan_up_6$Gene), length(expan_down_6$Gene),
                         length(rst_up_6$Gene), length(rst_down_6$Gene),
                         length(expan_up_12$Gene), length(expan_down_12$Gene),
                         length(rst_up_12$Gene), length(rst_down_12$Gene),
                         length(expan_up_24$Gene), length(expan_down_24$Gene),
                         length(rst_up_24$Gene), length(rst_down_24$Gene),
                         length(expan_up_48$Gene), length(expan_down_48$Gene),
                         length(rst_up_48$Gene), length(rst_down_48$Gene)))

de_norm_rep$DE_type = as.character(de_norm_rep$DE_type)
de_norm_rep$Reg_type = sapply(de_norm_rep$DE_type, function(x) {strsplit(x, ' ')[[1]][1]})
de_norm_rep$Reg_type[de_norm_rep$Reg_type == 'Down'] = 'Downregulated genes'
de_norm_rep$Reg_type[de_norm_rep$Reg_type == 'Up'] = 'Upregulated genes'

deg = ggplot(de_norm_rep, aes(reorder(DE_type, as.numeric(rownames(de_norm_rep))), Repeats_num_near_deg*100,
                              label = Number_of_repeats)) +
  geom_rect(data = de_norm_rep,
            aes(fill = Reg_type),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.2) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', fill = brewer.pal(8, 'Pastel1')[7]) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 40, hjust = 1),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        strip.text.x = element_text(size = 18,colour = "black"),
        strip.text.y = element_text(size = 18,colour = "black"),
        legend.position = 'none'
       # legend.text = element_text(size=18,color = "black"),
      #  legend.title = element_text(size = 18,color = "black")
      ) +
  labs(title = 'Amount of LINE/CR1 nearby differentially expressed regeneration genes',
       x = '', y = 'Percentage of repeats per one gene',
       fill = 'Repeat type') +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2])) +
  facet_grid(cols = vars(de_norm_rep$Repeat_type), rows = vars(Reg_type)) +
  geom_text(size = 5, hjust = .5, vjust = 0)


#ggsave(plot = deg, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/deg_all_regeneration_stages.tiff',
 #      width = 10, height = 10, device = 'tiff')


#----------table with ratios of CR1 repeats------------

de_norm_rep_up = subset(de_norm_rep, Reg_type == 'Upregulated genes')
str(de_norm_rep_up)
de_norm_rep_up = as.data.frame(t(de_norm_rep_up[, 1:3]))
colnames(de_norm_rep_up) = c("Up 0.5h exp", "Up 0.5h rest", "Up 3h exp", "Up 3h rest",
                             "Up 6h exp", "Up 6h rest", "Up 12h exp", "Up 12h rest",
                             "Up 24h exp", "Up 24h rest", "Up 48h exp", "Up 48h rest")
de_norm_rep_up = de_norm_rep_up[-1,]
de_norm_rep_up = de_norm_rep_up[-2,]
de_norm_rep_up = as.data.frame(t(sapply(de_norm_rep_up, function(x) as.numeric(x))))

ratio_up = data.frame(Ratio = c(de_norm_rep_up$`Up 0.5h rest`/de_norm_rep_up$`Up 0.5h exp`,
                      de_norm_rep_up$`Up 3h rest`/de_norm_rep_up$`Up 3h exp`,
                      de_norm_rep_up$`Up 6h rest`/de_norm_rep_up$`Up 6h exp`,
                      de_norm_rep_up$`Up 12h rest`/de_norm_rep_up$`Up 12h exp`,
                      de_norm_rep_up$`Up 24h rest`/de_norm_rep_up$`Up 24h exp`,
                      de_norm_rep_up$`Up 48h rest`/de_norm_rep_up$`Up 48h exp`))
rownames(ratio_up) = c('Up 0.5h', 'Up 3h', 'Up 6h', 'Up 12h', 'Up 24h', 'Up 48h')

de_norm_rep_down = subset(de_norm_rep, Reg_type == 'Downregulated genes')
str(de_norm_rep_down)
de_norm_rep_down = as.data.frame(t(de_norm_rep_down[, 1:3]))
colnames(de_norm_rep_down) = c("Down 0.5h exp", "Down 0.5h rest", "Down 3h exp", "Down 3h rest",
                             "Down 6h exp", "Down 6h rest", "Down 12h exp", "Down 12h rest",
                             "Down 24h exp", "Down 24h rest", "Down 48h exp", "Down 48h rest")
de_norm_rep_down = de_norm_rep_down[-1,]
de_norm_rep_down = de_norm_rep_down[-2,]
de_norm_rep_down = as.data.frame(t(sapply(de_norm_rep_down, function(x) as.numeric(x))))

ratio_down = data.frame(Ratio = c(de_norm_rep_down$`Down 0.5h rest`/de_norm_rep_down$`Down 0.5h exp`,
                                   de_norm_rep_down$`Down 3h rest`/de_norm_rep_down$`Down 3h exp`,
                                   de_norm_rep_down$`Down 6h rest`/de_norm_rep_down$`Down 6h exp`,
                                   de_norm_rep_down$`Down 12h rest`/de_norm_rep_down$`Down 12h exp`,
                                   de_norm_rep_down$`Down 24h rest`/de_norm_rep_down$`Down 24h exp`,
                                   de_norm_rep_down$`Down 48h rest`/de_norm_rep_down$`Down 48h exp`))
rownames(ratio_down) = c('Down 0.5h', 'Down 3h', 'Down 6h', 'Down 12h', 'Down 24h', 'Down 48h')

ratio_up_down = rbind(ratio_up, ratio_down)

ratio_up_down$Status = sapply(rownames(ratio_up_down), function(x) {strsplit(x, ' ')[[1]][1]})
ratio_up_down$Time = sapply(rownames(ratio_up_down), function(x) {strsplit(x, ' ')[[1]][2]})

ratio_up_down$Order = c(1:12)

str(ratio_up_down)

ratio_plot = ggplot(ratio_up_down, aes(x = reorder(Time, Order), y = Ratio, color = Status)) +
  geom_point(size = 4) +
  geom_line(aes(group = Status)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 18, face = 'bold'),
        axis.text.x = element_text(size=18,color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        legend.position = 'right',
        legend.text = element_text(size=18,color = "black"),
        legend.title = element_text(size=18,color = "black")
        # legend.text = element_text(size=18,color = "black"),
        #  legend.title = element_text(size = 18,color = "black")
  ) +
  labs(title = 'Rest/Expansion LINE/CR1 ratio nearby DEG',
       x = 'Time', y = 'Ratio',
       color = 'Expression') +
  ylim(c(-1,3)) +
  scale_color_manual(values = c(brewer.pal(11, 'RdBu')[9], brewer.pal(11, 'RdBu')[2]))

#ggsave(plot = ratio_plot,
 #      filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/ratio_plot.tiff',
  #     width = 10, height = 10, device = 'tiff')


#------------3 hours genes--------------------
length(unique(rst_down_3$Gene))
rst_down_3_genes = unique(rst_down_3$Gene)
write.csv(rst_down_3_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/rst_down_3_genes_CR1.csv")

length(unique(rst_up_3$Gene))
rst_up_3_genes = unique(rst_up_3$Gene)
write.csv(rst_up_3_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/rst_up_3_genes_CR1.csv")

length(unique(expan_down_3$Gene))
expan_down_3_genes = unique(expan_down_3$Gene)
write.csv(expan_down_3_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/expan_down_3_genes_CR1.csv")


length(unique(expan_up_3$Gene))
expan_up_3_genes = unique(expan_up_3$Gene)
write.csv(expan_up_3_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/expan_up_3_genes_CR1.csv")


#----------6 hours genes------------------------
length(unique(expan_up_6$Gene))
expan_up_6_genes = unique(expan_up_6$Gene)
write.csv(expan_up_6_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/expan_up_6_genes_CR1.csv")


length(unique(rst_up_6$Gene))
rst_up_6_genes = unique(rst_up_6$Gene)
write.csv(rst_up_6_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/rst_up_6_genes_CR1.csv")

#----------0.5 hours genes------------------------
length(unique(expan_up_05$Gene))
expan_up_05_genes = unique(expan_up_05$Gene)
write.csv(expan_up_05_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/expan_up_05_genes_CR1.csv")


length(unique(rst_up_05$Gene))
rst_up_05_genes = unique(rst_up_05$Gene)
write.csv(rst_up_05_genes,
          file = "/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/diff_expression_regeneration/rst_up_05_genes_CR1.csv")
#-------------GO analysis-------------------------------------

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
GO_3_up_expan$Term[GO_3_up_expan$Term=='determination of dorsal/ventral asymmetr...'] = 'determination of dorsal/ventral asymmetry'
GO_3_up_expan$Term[GO_3_up_expan$Term=='mammillothalamic axonal tract developmen...'] = 'mammillothalamic axonal tract development'

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

#GO 3h up rest

GO_3_up_rest = GO_analysis(rst_up_3)

plot_go_up_rest_3 = ggplot(GO_3_up_rest, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher),
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
  labs(title = 'Gene ontology of upregulated (3h) genes nearby non-expanded repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(GO_3_up_rest$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4))

ggsave(plot = plot_go_up_rest_3, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/plot_go_up_rest_3.tiff',
       width = 10, height = 10, device = 'tiff')  

#----------------6 hours stage--------------------------------

# GO 6h expan up
GO_6_up_expan = GO_analysis(expan_up_6)

GO_6_up_expan$Term[GO_6_up_expan$Term == 'signaling receptor ligand precursor proc...'] = 'signaling receptor ligand precursor processing'
GO_6_up_expan$Term[GO_6_up_expan$Term == 'positive regulation of cellular metaboli...'] = 'positive regulation of cellular metabolic process'

plot_go_up_expan_6 = ggplot(GO_6_up_expan, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher),
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
  labs(title = 'Gene ontology of upregulated (6h) genes nearby expanded repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(GO_6_up_expan$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4))

ggsave(plot = plot_go_up_expan_6, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/plot_go_up_expan_6.tiff',
       width = 10, height = 10, device = 'tiff')  

#GO 6h up rest

GO_6_up_rest = GO_analysis(rst_up_6)

GO_6_up_rest$Term[GO_6_up_rest$Term == 'signaling receptor ligand precursor proc...'] = 'signaling receptor ligand precursor processing'
GO_6_up_rest$Term[GO_6_up_rest$Term == 'regulation of angiotensin levels in bloo...'] = 'regulation of angiotensin levels in blood'
GO_6_up_rest$Term[GO_6_up_rest$Term == 'cellular response to environmental stimu...'] = 'cellular response to environmental stimulus'

plot_go_up_rest_6 = ggplot(GO_6_up_rest, aes(reorder(Term, -log10(classicFisher)), -log10(classicFisher),
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
  labs(title = 'Gene ontology of upregulated (6h) genes nearby non-expanded repeats',
       x = '', y = '-log10(p-value)',
       fill = 'Gene number') +
  scale_fill_gradient2(low = brewer.pal(9, name = "Set1")[2],
                       mid = 'white',
                       high = brewer.pal(9, name = "Set1")[1],
                       midpoint = mean(GO_6_up_rest$Significant)) +
  geom_text(aes(label = Significant),
            size = 6, hjust=1, vjust=1) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4))

ggsave(plot = plot_go_up_rest_6, filename = '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/plot_go_up_rest_6.tiff',
       width = 10, height = 10, device = 'tiff')  
