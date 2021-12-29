library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(gg.gap)
library(ggtree)
library(ape)
library(tidyverse)

#--------------SLR Hydra vulgaris----------------------

df_brown = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/ete3_SLR/SLR~0c474aaf0a64db71469f60856e944290/out'
                       )

col_names = gsub("Adj.Pval", "Adj_Pval", colnames(df_brown)) %>%
  gsub("X", "", .) %>%
  gsub("Q.value", "Q_value", .)

col_names = unlist(strsplit(col_names, '[.]'))
col_names = unlist(lapply(col_names, function(x){x[!x ==""]}))

df_brown = gsub("^   ","",df_brown$X..Site..Neutral..Optimal...Omega....lower....upper.LRT_Stat....Pval.....Adj.Pval....Q.value.Result.Note) %>% 
  gsub("^    ","",.) %>% 
  gsub("^   ","",.) %>% 
  gsub("^  ","",.) %>%
  gsub("^ ","",.) %>%
  gsub("     "," ",.) %>%
  gsub("    "," ",.) %>%
  gsub("   "," ",.) %>% 
  gsub("!","",.) %>%
  gsub("  "," ",.) %>%
  gsub("Single char","Single_char",.) %>%
  gsub("All gaps","All_gaps",.) %>%
 # gsub("\\(.*\\)","",.)%>%
  as.data.frame()

df_brown <- stringr::str_split_fixed(df_brown$.," ",Inf) %>% as.data.frame()

colnames(df_brown) = col_names

str(df_brown)
df_brown$Omega = as.numeric(df_brown$Omega)
df_brown$Adj_Pval = as.numeric(df_brown$Adj_Pval)

df_brown$Significant = apply(df_brown, 1,
                             function(x){
                            ifelse(x[11] == '' , 'notsign',
                            ifelse(x[11] == '++', 'positive', 'negative'))
                              }
                               )
df_brown$Significant = as.factor(df_brown$Significant)
df_brown$Omega[df_brown$Significant == 'notsign'] = 0

plot = ggplot(df_brown, aes(as.numeric(rownames(df_brown)), Omega, fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray9') +
  theme_bw() +
  geom_hline(yintercept = c(1),linetype="dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
        ) +
  labs(title = 'dN/dS values for Hydra vulgaris CR1 reverse transcriptase sites',
       x = 'Position', y = 'dN/dS'
  ) +
 # scale_y_continuous(breaks=seq(0, 70, 1), limits=c(0, 75)) +
  scale_x_continuous(breaks=seq(1, 284), limits=c(0, 285)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9],'white', brewer.pal(11, 'RdBu')[2]))

plot2 = gg.gap(plot = plot,
       ylim = c(0,75),
       segments = c(2, 70))

ggsave(plot = plot2, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvul_SLR.tiff',
       device = 'tiff', width = 10, height = 10)

#-------------SLR Hydra viridissima-------------------------------
df_green = read.delim2('/Users/alenakizenko/Documents/PhD/project/green_hydra_repeat_ann/rvt1_analysis/ete3_SLR/SLR~3af5025cd4252949b07c2abd461d9296/out'
)

col_names = gsub("Adj.Pval", "Adj_Pval", colnames(df_green)) %>%
  gsub("X", "", .) %>%
  gsub("Q.value", "Q_value", .)

col_names = unlist(strsplit(col_names, '[.]'))
col_names = unlist(lapply(col_names, function(x){x[!x ==""]}))

df_green = gsub("^   ","",df_green$X..Site..Neutral..Optimal...Omega....lower....upper.LRT_Stat....Pval.....Adj.Pval....Q.value.Result.Note) %>% 
  gsub("^    ","",.) %>% 
  gsub("^   ","",.) %>% 
  gsub("^  ","",.) %>%
  gsub("^ ","",.) %>%
  gsub("     "," ",.) %>%
  gsub("    "," ",.) %>%
  gsub("   "," ",.) %>% 
  gsub("!","",.) %>%
  gsub("  "," ",.) %>%
  gsub("Single char","Single_char",.) %>%
  gsub("All gaps","All_gaps",.) %>%
  # gsub("\\(.*\\)","",.)%>%
  as.data.frame()

df_green <- stringr::str_split_fixed(df_green$.," ",Inf) %>% as.data.frame()

colnames(df_green) = col_names

str(df_green)
df_green$Omega = as.numeric(df_green$Omega)
df_green$Adj_Pval = as.numeric(df_green$Adj_Pval)
df_green$Result[df_green$Result == 'Constant'] = ''

df_green$Significant = apply(df_green, 1,
                             function(x){
                               ifelse(x[11] == '' , 'notsign',
                                      ifelse(x[11] == '++', 'positive', 'negative'))
                             }
)
df_green$Significant = as.factor(df_green$Significant)
df_green$Omega[df_green$Significant == 'notsign'] = 0

plot = ggplot(df_green, aes(as.numeric(rownames(df_green)), Omega, fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray9') +
  theme_classic() +
  geom_hline(yintercept = c(1),linetype="dashed") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
       ) +
  labs(title = 'dN/dS values for Hydra viridissima CR1 reverse transcriptase sites',
       x = 'Position', y = 'dN/dS'
  ) +
   scale_y_continuous(breaks=seq(0, 3, 1), limits=c(0, 4)) +
  scale_x_continuous(breaks=seq(1, 234), limits=c(0, 235)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], 'white', brewer.pal(11, 'RdBu')[2]))


ggsave(plot = plot, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvir_SLR.tiff',
       device = 'tiff', width = 10, height = 10)

#-------------SLR Hydra oligactis-------------------------------
df_holi = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_oligactis/rt_fr/SLR~c64180e266bfcf1f1b619983cf016c91/out'
)

col_names = gsub("Adj.Pval", "Adj_Pval", colnames(df_holi)) %>%
  gsub("X", "", .) %>%
  gsub("Q.value", "Q_value", .)

col_names = unlist(strsplit(col_names, '[.]'))
col_names = unlist(lapply(col_names, function(x){x[!x ==""]}))

df_holi = gsub("^   ","",df_holi$X..Site..Neutral..Optimal...Omega....lower....upper.LRT_Stat....Pval.....Adj.Pval....Q.value.Result.Note) %>% 
  gsub("^    ","",.) %>% 
  gsub("^   ","",.) %>% 
  gsub("^  ","",.) %>%
  gsub("^ ","",.) %>%
  gsub("!","",.) %>%
 # gsub("     "," ",.) %>%
  gsub("    "," ",.) %>%
  gsub("   "," ",.) %>%
  gsub("  "," ",.) %>%
  gsub("Single char","",.) %>%
  gsub("All gaps","",.) %>%
  gsub("Synonymous","",.) %>%
  gsub("Constant","",.) %>%
  # gsub("\\(.*\\)","",.)%>%
  as.data.frame()

df_holi <- stringr::str_split_fixed(df_holi$.," ",Inf) %>% as.data.frame()

colnames(df_holi) = col_names

str(df_holi)
df_holi$Omega = as.numeric(df_holi$Omega)
df_holi$Adj_Pval = as.numeric(df_holi$Adj_Pval)

df_holi$Significant = apply(df_holi, 1,
                             function(x){
                               ifelse(x[11] == '' , 'notsign', 'negative')
                             }
)
df_holi$Significant = as.factor(df_holi$Significant)
df_holi$Omega[df_holi$Significant == 'notsign'] = 0

plot = ggplot(df_holi, aes(as.numeric(rownames(df_holi)), Omega, fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray9') +
  theme_classic() +
  geom_hline(yintercept = c(1),linetype="dashed") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title = 'dN/dS values for Hydra oligactis CR1 reverse transcriptase sites',
       x = 'Position', y = 'dN/dS'
  ) +
  scale_y_continuous(breaks=seq(0, 3, 1), limits=c(0, 4)) +
  scale_x_continuous(breaks=seq(1, 283), limits=c(0, 284)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], 'white', brewer.pal(11, 'RdBu')[2]))


ggsave(plot = plot, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_holi_SLR.tiff',
       device = 'tiff', width = 10, height = 10)


#-------------fb model Hydra vulgaris-------------------------------


fb_tree_hvul = read.tree('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/ete3_fb/out_tree_prefix.fasta')

tip_labels = fb_tree_hvul$tip.label
fb_tree_hvul = root(fb_tree_hvul, 'old_hvul_LINE/CR1_rnd-1_family-110_117')
grp <- list(
  'Hydra vulgaris (new)' = tip_labels[grepl("new", tip_labels)],
  'Hydra vulgaris (old)' = tip_labels[grepl("old", tip_labels)])


tree2 = groupOTU(fb_tree_hvul, grp)

tree2$edge.length = ifelse(tree2$edge.length <= 1, 'positive', 'negative')

p = ggtree(tree2,
         #  layout = 'circular',
           aes(color=branch.length),
           branch.length='none'
) +
#  geom_text(aes(label=branch.length), size = 3) +
  geom_tippoint(aes(shape = group, color= group), size=2, alpha = .8) +
  #scale_color_gradient2(low = brewer.pal(9, name = "Set1")[2],
   #                   # mid = 'white',
    #                   high = brewer.pal(9, name = "Set1")[1],
     #                  midpoint =200,
      #                 ) +
  scale_color_manual(values = c('white',
                                brewer.pal(9, name = "BrBG")[1],
                                brewer.pal(9, name = "BrBG")[3],
                                brewer.pal(9, name = "Set1")[2],
                                brewer.pal(9, name = "Set1")[1])) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_fb.tiff',
       width = 10, height = 10, device = 'tiff')


#-------------fb model Hydra viridissima-------------------------------


fb_tree_hvir = read.tree('/Users/alenakizenko/Documents/PhD/project/green_hydra_repeat_ann/rvt1_analysis/ete3_fb/fb~addefe208719013c6c83b4e8ddb7e6e8/out_tree')

tip_labels = fb_tree_hvir$tip.label
fb_tree_hvir = root(fb_tree_hvir, 'hvir_LINE/CR1_rnd-1_family-457_56')

fb_tree_hvir$edge.length = ifelse(fb_tree_hvir$edge.length <= 1, 'positive', 'negative')

p = ggtree(fb_tree_hvir,
           #  layout = 'circular',
           aes(color=branch.length),
           branch.length='none'
) +
  #geom_text(aes(label=branch.length), size = 3) +
  geom_tippoint(color = brewer.pal(3, 'Dark2')[1],size=2, alpha = .8) +
  #scale_color_gradient2(low = brewer.pal(9, name = "Set1")[2],
  #                   # mid = 'white',
  #                   high = brewer.pal(9, name = "Set1")[1],
  #                  midpoint =200,
  #                 ) +
  scale_color_manual(values = c('white',
                                brewer.pal(9, name = "Set1")[2],
                                brewer.pal(9, name = "Set1")[1])) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_fb_hvir.tiff',
       width = 10, height = 10, device = 'tiff')


#-----b_neut and b_free statistics hvir-----------------

out_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/green_hydra_repeat_ann/rvt1_analysis/ete3_b_free_b_neut/b_free_b_neut_out.tsv',
                       sep = '\t', header = F)

colnames(out_hvir) = c('Model', 'Branch', 'lnL', 'np')
  
b_neut = subset(out_hvir, Model == 'b_neut')
b_free = subset(out_hvir, Model == 'b_free')
b_neut_b_free = merge(b_neut, b_free, by = 'Branch')
str(b_neut_b_free)
options(digits=10)
b_neut_b_free$lnL.x = as.numeric(b_neut_b_free$lnL.x)
b_neut_b_free$lnL.y = as.numeric(b_neut_b_free$lnL.y)
b_neut_b_free$np.x = as.numeric(b_neut_b_free$np.x)
b_neut_b_free$np.y = as.numeric(b_neut_b_free$np.y)

b_neut_b_free$p_value = 1 - pchisq(2*(b_neut_b_free$lnL.y - b_neut_b_free$lnL.x),
                               df=(b_neut_b_free$np.y - b_neut_b_free$np.x))

#-----final tree hvir-------


bfree_tree_hvir = read.tree('/Users/alenakizenko/Documents/PhD/project/hydra_viridissima/rvt1_analysis/ete3_b_free_b_neut/output_tree_sign')

tip_labels = bfree_tree_hvir$tip.label
bfree_tree_hvir$edge.length = ifelse(bfree_tree_hvir$edge.length > 1, 'directional',
                          ifelse(bfree_tree_hvir$edge.length == 0, 'purifying', 'neutral'))

p = ggtree(bfree_tree_hvir,
           layout = 'circular',
           color='black',
           branch.length='none'
) +
  #  geom_text(aes(label=branch.length), size = 3) +
  geom_tippoint(aes(color = as.factor(branch.length)), size=4, alpha = .8) +
  scale_color_manual(values = c(brewer.pal(11, 'RdBu')[9], 'gray55', brewer.pal(11, 'RdBu')[2]),
                     labels = c('Purifying', 'Neutral', 'Directional')) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(shape="Repeat type",col="Selection")
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_b_free_b_neut.tiff',
       width = 10, height = 10, device = 'tiff')
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_bfree_hvir.tiff',
       width = 10, height = 10, device = 'tiff')

#------------transposase DNA/TcMar-Tc1--------------

df_transp_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/transposase_omega/out_hvir'
)

col_names = gsub("Adj.Pval", "Adj_Pval", colnames(df_transp_hvir)) %>%
  gsub("X", "", .) %>%
  gsub("Q.value", "Q_value", .)

col_names = unlist(strsplit(col_names, '[.]'))
col_names = unlist(lapply(col_names, function(x){x[!x ==""]}))

df_transp_hvir = gsub("^   ","",df_transp_hvir$X..Site..Neutral..Optimal...Omega....lower....upper.LRT_Stat....Pval.....Adj.Pval....Q.value.Result.Note) %>% 
  gsub("^    ","",.) %>% 
  gsub("^   ","",.) %>% 
  gsub("^  ","",.) %>%
  gsub("^ ","",.) %>%
  gsub("!","",.) %>%
  # gsub("     "," ",.) %>%
  gsub("    "," ",.) %>%
  gsub("   "," ",.) %>%
  gsub("  "," ",.) %>%
  gsub("Single char","",.) %>%
  gsub("All gaps","",.) %>%
  gsub("Synonymous","",.) %>%
  gsub("Constant","",.) %>%
  # gsub("\\(.*\\)","",.)%>%
  as.data.frame()

df_transp_hvir <- stringr::str_split_fixed(df_transp_hvir$.," ",Inf) %>% as.data.frame()

colnames(df_transp_hvir) = col_names

str(df_transp_hvir)
df_transp_hvir$Omega = as.numeric(df_transp_hvir$Omega)
df_transp_hvir$Adj_Pval = as.numeric(df_transp_hvir$Adj_Pval)

df_transp_hvir$Significant = apply(df_transp_hvir, 1,
                            function(x){
                              ifelse(x[11] == '' , 'notsign', ifelse(x[11] == '+', 'positive', 'negative'))
                            }
)
df_transp_hvir$Significant = as.factor(df_transp_hvir$Significant)
df_transp_hvir$Omega[df_transp_hvir$Significant == 'notsign'] = 0

plot = ggplot(df_transp_hvir, aes(as.numeric(rownames(df_transp_hvir)), Omega, fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray9') +
  theme_bw() +
  geom_hline(yintercept = c(1),linetype="dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title = 'dN/dS values for Hydra vulgaris DNA/TcMar-Tc1 reverse transcriptase sites',
       x = 'Position', y = 'dN/dS'
  ) +
  # scale_y_continuous(breaks=seq(0, 70, 1), limits=c(0, 75)) +
  scale_x_continuous(breaks=seq(1, 54), limits=c(0, 55)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9],'white', brewer.pal(11, 'RdBu')[2]))

plot2 = gg.gap(plot = plot,
               ylim = c(0,100),
               segments = c(2, 98))

ggsave(plot = plot2, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvir_SLR_DNA.tiff',
       device = 'tiff', width = 10, height = 10)

#-----------hvul DNA/TcMar-Tc1---------------------------

df_transp_hvir = read.delim2('/Users/alenakizenko/Documents/PhD/project/transposase_omega/out_hvul'
)

col_names = gsub("Adj.Pval", "Adj_Pval", colnames(df_transp_hvir)) %>%
  gsub("X", "", .) %>%
  gsub("Q.value", "Q_value", .)

col_names = unlist(strsplit(col_names, '[.]'))
col_names = unlist(lapply(col_names, function(x){x[!x ==""]}))

df_transp_hvir = gsub("^   ","",df_transp_hvir$X..Site..Neutral..Optimal...Omega....lower....upper.LRT_Stat....Pval.....Adj.Pval....Q.value.Result.Note) %>% 
  gsub("^    ","",.) %>% 
  gsub("^   ","",.) %>% 
  gsub("^  ","",.) %>%
  gsub("^ ","",.) %>%
  gsub("!","",.) %>%
  # gsub("     "," ",.) %>%
  gsub("    "," ",.) %>%
  gsub("   "," ",.) %>%
  gsub("  "," ",.) %>%
  gsub("Single char","",.) %>%
  gsub("All gaps","",.) %>%
  gsub("Synonymous","",.) %>%
  gsub("Constant","",.) %>%
  # gsub("\\(.*\\)","",.)%>%
  as.data.frame()

df_transp_hvir <- stringr::str_split_fixed(df_transp_hvir$.," ",Inf) %>% as.data.frame()

colnames(df_transp_hvir) = col_names

str(df_transp_hvir)
df_transp_hvir$Omega = as.numeric(df_transp_hvir$Omega)
df_transp_hvir$Adj_Pval = as.numeric(df_transp_hvir$Adj_Pval)

df_transp_hvir$Significant = apply(df_transp_hvir, 1,
                                   function(x){
                                     ifelse(x[11] == '' , 'notsign', ifelse(x[11] == '+', 'positive', 'negative'))
                                   }
)
df_transp_hvir$Significant = as.factor(df_transp_hvir$Significant)
df_transp_hvir$Omega[df_transp_hvir$Significant == 'notsign'] = 0

plot = ggplot(df_transp_hvir, aes(as.numeric(rownames(df_transp_hvir)), Omega, fill = Significant)) +
  geom_bar(stat = 'identity', color = 'gray9') +
  theme_classic() +
  geom_hline(yintercept = c(1),linetype="dashed") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=18,color = "black"),
        axis.title.x = element_text(size=18,color = "black"),
        axis.title.y = element_text(size=18,color = "black"),
        #strip.text.x = element_text(size = 18,colour = "black"),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title = 'dN/dS values for Hydra viridissima DNA/TcMar-Tc1 transposase sites',
       x = 'Position', y = 'dN/dS'
  ) +
  scale_y_continuous(breaks=seq(0, 3, 1), limits=c(0, 4)) +
  scale_x_continuous(breaks=seq(1, 62), limits=c(0, 63)) +
  scale_fill_manual(values = c(brewer.pal(11, 'RdBu')[9], 'white', brewer.pal(11, 'RdBu')[2]))


ggsave(plot = plot, filename = '/Users/alenakizenko/Documents/PhD/project/pictures/dN_dS_hvul_SLR_DNA.tiff',
       device = 'tiff', width = 10, height = 10)

#------------b_neut_b_free hvul-----------

b_free_b_neut_tree_hvul = read.tree('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/ete3_b_neut_b_free/tree_prefix.fasta')

tip_labels = b_free_b_neut_tree_hvul$tip.label

#b_free_b_neut_tree_hvul = root(b_free_b_neut_tree_hvul, 'old_hvul_LINE/CR1_rnd-1_family-110_117')
grp <- list(
  'Bubble' = tip_labels[grepl("new", tip_labels)],
  'Half moon' = tip_labels[grepl("old", tip_labels)])


tree2 = groupOTU(b_free_b_neut_tree_hvul, grp)

#tree2$edge.length = ifelse(tree2$edge.length > 1, 'directional',
 #                          ifelse(tree2$edge.length == 0, 'purifying', 'neutral'))

p = ggtree(tree2,
           layout = 'circular',
           color='black',
           branch.length='none'
) +
  #  geom_text(aes(label=branch.length), size = 3) +
  geom_tippoint(aes(shape = group, color = as.factor(branch.length)), size=4, alpha = .8) +
  scale_color_manual(values = c(brewer.pal(11, 'RdBu')[9], 'gray55', brewer.pal(11, 'RdBu')[2]),
                    labels = c('Purifying', 'Neutral', 'Directional')) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(shape="Repeat type",col="Selection")
p

ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_b_free_b_neut.tiff',
       width = 10, height = 10, device = 'tiff')


#------------tree with labesl-----------

lab_tree = read.tree('/Users/alenakizenko/Documents/PhD/project/rvt1_hvir_hvul/hvul_hvir_rvt1_CR1_codon_align_final_lab.treefile')
lab_tree = root(lab_tree, "hvir_LINE_CR1_rnd-1_family-457_43_-_")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols = sample(col_vector, 10)

p = ggtree(lab_tree,
           layout = 'circular',
           color='gray66'
           #branch.length='none'
) +
  #  geom_text(aes(label=branch.length), size = 3) +
  #geom_tippoint(size=4, alpha = .8) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  geom_hilight(node=which(lab_tree$node.label == '#0'),
               fill=cols[1], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#1'),
               fill=cols[2], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#2'),
               fill=cols[3], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#3'),
               fill=cols[4], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#4'),
               fill=cols[5], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#5'),
               fill=cols[6], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#6'),
               fill=cols[7], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#7'),
               fill=cols[8], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#8'),
               fill=cols[9], alpha=1, extend = .5) +
  geom_hilight(node=which(lab_tree$node.label == '#10'),
               fill=cols[10], alpha=1, extend = .5) +
  geom_nodelab(size = 2)
p

# ----------------b_neut b_free internal-----------

tree <- read.tree('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/RT1_analysis/ete3_b_neut_b_free_internal/result/tree_for_plot/hvul_hvir_rvt1_CR1_codon_align.fasta_prefix_branch_omega.fasta')
tip_labels = tree$tip.label
edge_length = tree$edge.length

grp <- list(
  'Hydra viridissima' = tip_labels[grepl("hvir", tip_labels)],
  'Hydra vulgaris new' = tip_labels[grepl("new", tip_labels)],
  'Hydra vulgaris old' = tip_labels[grepl("old", tip_labels)])


tree2 = groupOTU(tree, grp)

p = ggtree(tree2,
           layout = 'circular',
           aes(subset = !isTip, color=as.numeric(label)),
           size = 1
) +
  scale_colour_gradient(low = 'white', high = 'black') +
  #geom_tippoint(aes(shape = group),
   #             color = 'green',
    #            fill = 'white',
     #           size=4, alpha = .8) +
  theme(
    plot.title = element_blank(),
    legend.position="left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  #geom_nodelab(size = 2) +
  geom_hilight(node = 157, fill=brewer.pal(11, 'BrBG')[1], alpha=.6, extend = .5) +
  geom_hilight(node=238, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=105, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=138, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=257, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=240, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=252, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=106, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=125, fill=brewer.pal(11, 'BrBG')[1], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=126, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=127, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=1, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) +
  geom_hilight(node=144, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_balance(node=141, fill=brewer.pal(11, 'BrBG')[3], color='white', alpha=0.6, extend=.5) +
  geom_balance(node=143, fill=brewer.pal(11, 'BrBG')[10], color='white', alpha=0.6, extend=.5) 
p


ggsave(plot = p,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/RT_b_free_b_neut_internal.tiff',
       width = 10, height = 10, device = 'tiff')

