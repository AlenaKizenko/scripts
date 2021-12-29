library(ggplot2)
library(RColorBrewer)

cnidaria = data.frame(Species = c('Aurelia sp', 'Trachythela sp.', 'Cassiopea andromeda', 'Pachyseris speciosa', 'Anthopleura sola',
                                  'Chrysaora achlyos', 'Alatinidae', 'Stichodactyla helianthus', 'Acropora intermedia', 'Acropora selago',
                                  'Acropora microphthalma', 'Acropora gemmifera', 'Acropora echinata', 'Acropora cytherea', 'Acropora awi',
                                  'Acropora acuminata', 'Pocillopora verrucosa', 'Exaiptasia diaphana', 'Sanderia malayensis', 'Heteractis crispa',
                                  'Stichodactyla mertensii', 'Heteractis magnifica', 'Chrysaora chesapeakei', 'Aurelia coerulea', 'Myxobolus squamalis',
                                  'Carybdea marsupialis', 'Chrysaora fuscescens', 'Henneguya salminicola', 'Actinia tenebrosa', 'Montipora capitata',
                                  'Montipora efflorescens', 'Actinia equina', 'Hydra viridissima', 'Morbakka virulenta', 'Porites rus',
                                  'Calvadosia cruxmelitensis', 'Cassiopea xamachana', 'Nemopilema nomurai', 'Rhopilema esculentum', 'Renilla reniformis',
                                  'Enteromyxum leei', 'Sphaeromyxa zaharoni', 'Phymanthus crucifer', 'Kudoa iwatai', 'Astreopora myriophthalma',
                                  'Acropora muricata', 'Acropora yongei', 'Acropora nasuta', 'Acropora florida', 'Pocillopora damicornis',
                                  'Clytia hemisphaerica', 'Thelohanellus kitauei', 'Chrysaora quinquecirrha', 'Acropora hyacinthus', 'Alatina alata',
                                  'Anemonia viridis', 'Orbicella faveolata', 'Hydra vulgaris', 'Stylophora pistillata', 'Craspedacusta sowerbii',
                                  'Acropora digitifera', 'Dendronephthya gigantea', 'Hydra oligactis', 'Aurelia aurita', 'Montipora cactus',
                                  'Acropora tenuis', 'Acropora millepora', 'Nematostella vectensis'),
                      Taxon = c('Scyphozoa', 'Anthozoa', 'Scyphozoa', 'Anthozoa', 'Anthozoa',
                                'Scyphozoa', 'Cubozoa', 'Anthozoa', 'Anthozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Anthozoa', 'Anthozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Anthozoa', 'Scyphozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Scyphozoa','Scyphozoa', 'Myxozoa',
                                'Cubozoa', 'Scyphozoa', 'Myxozoa', 'Anthozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Hydrozoa', 'Cubozoa', 'Anthozoa',
                                'Staurozoa', 'Scyphozoa', 'Scyphozoa', 'Scyphozoa', 'Anthozoa',
                                'Myxozoa', 'Myxozoa', 'Anthozoa', 'Myxozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Anthozoa', 'Anthozoa', 'Anthozoa',
                                'Hydrozoa', 'Myxozoa', 'Scyphozoa', 'Anthozoa', 'Cubozoa',
                                'Anthozoa', 'Anthozoa', 'Hydrozoa', 'Anthozoa', 'Hydrozoa',
                                'Anthozoa', 'Anthozoa', 'Hydrozoa', 'Scyphozoa', 'Anthozoa',
                                'Anthozoa', 'Anthozoa', 'Anthozoa'),
                      Genome_size = c(426.02, 578.19, 406.64, 625.35, 9.81,
                                      432.41, 1400.7, 296.17, 416.88, 392.94,
                                      383.74, 400.99, 401.48, 426.32, 428.84,
                                      394.74, 380.5, 256.13, 184.37, 399.04,
                                      292.95, 376.95, 457.12, 669.42, 43.66,
                                      1286.92, 266.41, 61.44, 238.18, 614.5,
                                      643.32, 409.06, 284.27, 951.58, 469.75,
                                      209.38, 393.49, 213.62, 256.68, 131.56,
                                      68.1635, 74.58, 437.44, 31.18, 373.48,
                                      420.73, 438.05, 416.42, 442.83, 234.33,
                                      420.87, 150.35, 337.42, 447.2, 1075.36,
                                      400.6, 485.53, 1250, 400.1, 673.19,
                                      447.48, 286.13, 1450, 376.93, 652.74,
                                      403.14, 386.6, 356.61),
                      GC_content = c(32.1, 37.3, 33.7, 39.6, 33.8,
                                     42.1, 40.0, 37.9, 36.9, 36.5,
                                     35.3, 35.2, 33.2, 36.6, 33.8,
                                     36.9, 38.0, 38.8, 37.3, 39.1,
                                     35.8, 38.1, 40.2, 40.9, 27.3,
                                     38.8, 37.6, 29.0, 38.5, 36.8,
                                     36.0, 37.5, 24.7, 31.4, 0.0,
                                     0.0, 0.0, 38.4, 36.3, 36.6,
                                     48.686, 28.0, 43.3, 23.6, 38.4,
                                     36.4, 36.4, 36.2, 35.4, 38.3,
                                     35.1, 34.5, 35.6, 35.9, 38.0,
                                     38.3, 41.9, 29.1, 39.7, 40.8,
                                     40.5, 36.8, 28.0, 34.6, 36.4,
                                     36.1, 39.9, 41.9))

write.csv2(x = cnidaria, file = '/Users/alenakizenko/Documents/PhD/project/cnidaria_genomes.csv')

str(cnidaria)
cnidaria = subset(cnidaria, Taxon != 'Myxozoa')

cnidaria$Taxon = as.factor(cnidaria$Taxon)


cnidaria = cnidaria[order(cnidaria$Taxon),]
cnidaria$Order = c(1:length(cnidaria$Species))



cnid = ggplot(cnidaria, aes(x = reorder(Species, Order), y = Genome_size, fill = Taxon)) +
  geom_bar(stat = 'identity', alpha = .9) +
  coord_polar('x') +
  theme_bw() +
  theme(legend.position = c(0.01, 0.2),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 18),
        #axis.line=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 18),
        #panel.background=element_blank(),
        panel.border =element_blank(),
        #panel.grid= element_line(color = 'gray66', size = .1),
        #panel.grid = element_blank(),
        #panel.grid.minor=element_blank(),
        #plot.background=element_blank(),
        ) +
  scale_fill_brewer(type = 'seq', palette = 'Spectral') +
  ylab(label = 'Genome size, bp')

cnid

ggsave(plot = cnid, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/cnidaria_genomes.tiff',
      device = 'tiff', width = 10, height = 10)

