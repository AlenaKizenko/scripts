#load library
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)


#-------------------------HYLI-------------------------------------------------
#read in entire BAM file
bam <- scanBam("/Users/alenakizenko/Documents/PhD/project/pirna/pirna_aep_analysis/hyli/hyli_align_sorted.bam")

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
# expansion
#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#-------------------------------expansion---------------------------------------
#store the mapped positions on the plus and minus strands
exp_neg <- bam_df[bam_df$rname == 'LINE_CR1_expansion' &
                      apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
]

length(exp)

exp_pos <- bam_df[bam_df$rname == 'LINE_CR1_expansion' &
                      apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
]
length(rst)

#calculate the densities
exp_neg_density <- density(exp_neg)
exp_pos_density <- density(exp_pos)
#display the negative strand with negative values
exp_neg_density$y <- exp_neg_density$y * -1

exp_final = data.frame('Coordinate_pos' = exp_pos_density$x, 'Coordinate_neg' = exp_neg_density$x,
                         'Density_pos' = exp_pos_density$y, 'Density_neg' = exp_neg_density$y)

exp_pirna_plot = ggplot(exp_final, aes(Coordinate_pos, Density_pos)) +
  geom_line(stat = 'identity', color = brewer.pal(11, 'Spectral')[1]) +
  geom_line(aes(Coordinate_neg, Density_neg), stat = 'identity',
           color = brewer.pal(11, 'Spectral')[10]) +
  geom_hline(yintercept = 0, color = 'gray55') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  #ylim(-8e-07, 8e-07) +
  labs(x="Hydra vulgaris specific expansion", y="Density")

ggsave(plot = exp_pirna_plot,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/hvulg_purna_exp.tiff',
       width = 10, height = 10, device = 'tiff')

#----------------------------rest-----------------------------------------------
#store the mapped positions on the plus and minus strands
rest_neg <- bam_df[bam_df$rname == 'LINE_CR1_rest' &
                    apply(as.data.frame(bam_df$flag), 1, check_neg),
                  'pos'
]

rest_pos <- bam_df[bam_df$rname == 'LINE_CR1_rest' &
                    apply(as.data.frame(bam_df$flag), 1, check_pos),
                  'pos'
]

#calculate the densities
rest_neg_density <- density(rest_neg)
rest_pos_density <- density(rest_pos)
#display the negative strand with negative values
rest_neg_density$y <- rest_neg_density$y * -1

rest_final = data.frame('Coordinate_pos' = rest_pos_density$x, 'Coordinate_neg' = rest_neg_density$x,
                       'Density_pos' = rest_pos_density$y, 'Density_neg' = rest_neg_density$y)
rest_pirna_plot = ggplot(rest_final, aes(Coordinate_pos, Density_pos)) +
  geom_line(stat = 'identity', color = brewer.pal(11, 'Spectral')[1]) +
  geom_line(aes(Coordinate_neg, Density_neg), stat = 'identity',
            color = brewer.pal(11, 'Spectral')[10]) +
  geom_hline(yintercept = 0, color = 'gray55') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
 # ylim(-8e-07, 8e-07) +
  labs(x="Brown hydra specific expansion", y="Density")

ggsave(plot = rest_pirna_plot,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/hvulg_purna_rest.tiff',
       width = 10, height = 10, device = 'tiff')

hist(rest_neg)
hist(rest_pos)
hist(exp_neg)
hist(exp_pos)
plot(density(exp_pos))

# idxstats hyli

# dividing sequences by bins

library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(dplyr)

lmig_hvul = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/rma1_replanscape/hvul_divsum.txt',
                       header = T, sep = '\t', skip = 4)
lmig_hvul = lmig_hvul[-1, ]
stop_num = which(lmig_hvul$Class == 'Coverage for each repeat class and divergence (Kimura)')
lmig_hvul = as.data.frame(lmig_hvul[c(1:stop_num-1),])
cr1_hvul = subset(lmig_hvul, Class == 'LINE/CR1' & Kimura. > 0)
str(cr1_hvul)
cr1_hvul$absLen = as.numeric(cr1_hvul$absLen)
cr1_hvul$wellCharLen = as.numeric(cr1_hvul$wellCharLen)
cr1_hvul$Kimura. = as.numeric(cr1_hvul$Kimura.)

ggplot(cr1_hvul, aes(Kimura.)) +
  geom_histogram(binwidth = 1, color = 'gray5', fill = 'gray55') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(6, 33), breaks = c(6:33))

tags <- c("[6-11)","[11-17)", "[17-23)", "[23-33)")
cr1_fam <- cr1_hvul %>% select(Repeat, Kimura.) #pick the variable 
cr1_binned <- as_tibble(cr1_fam) %>% 
  mutate(tag = case_when(
    Kimura. >= 6 & Kimura. < 11 ~ tags[1],
    Kimura. >= 11 & Kimura. < 17 ~ tags[2],
    Kimura. >= 17 & Kimura. < 23 ~ tags[3],
    Kimura. >= 23 & Kimura. < 33 ~ tags[4]
  ))

write(x = subset(cr1_binned, tag == '[6-11)')$Repeat, file = '/Users/alenakizenko/Documents/PhD/project/pirna/Kimura_binning/cr1_6_11.txt')
write(x = subset(cr1_binned, tag == '[11-17)')$Repeat, file = '/Users/alenakizenko/Documents/PhD/project/pirna/Kimura_binning/cr1_11_17.txt')
write(x = subset(cr1_binned, tag == '[17-23)')$Repeat, file = '/Users/alenakizenko/Documents/PhD/project/pirna/Kimura_binning/cr1_17_23.txt')
write(x = subset(cr1_binned, tag == '[23-33)')$Repeat, file = '/Users/alenakizenko/Documents/PhD/project/pirna/Kimura_binning/cr1_23_33.txt')

#-----------------------

df_hyli = read.delim('/Users/alenakizenko/Documents/PhD/project/pirna/pirna_aep_analysis/hyli/hyli_idxstats.tsv',
                sep = '\t', header = F)
df_hywi = read.delim('/Users/alenakizenko/Documents/PhD/project/pirna/pirna_aep_analysis/hywi/hywi_idxstats.tsv',
                     sep = '\t', header = F)

names(df_hyli) = c('Sequence_name', 'Sequence_length', 'Mapped_read-segments', 'Unmapped_read-segments')
df_hyli$Normalized_mapped = df_hyli$`Mapped_read-segments`/df_hyli$Sequence_length

names(df_hywi) = c('Sequence_name', 'Sequence_length', 'Mapped_read-segments', 'Unmapped_read-segments')
df_hywi$Normalized_mapped = df_hywi$`Mapped_read-segments`/df_hywi$Sequence_length

hywi_hyli_df = data.frame(Normalized_reads = c(df_hyli[-3,]$Normalized_mapped, df_hywi[-3,]$Normalized_mapped),
                          Repeat_type = c('New', 'Old', 'New', 'Old'),
                          Species = c('Hyli', 'Hyli', 'Hywi', 'Hywi'))

str(hywi_hyli_df)

hywi_hyli_df$Repeat_type = as.factor(hywi_hyli_df$Repeat_type)
hywi_hyli_df$Species = as.factor(hywi_hyli_df$Species)

read_plot_hyli_hywi = ggplot(hywi_hyli_df, aes(Species, Normalized_reads*100, fill = Repeat_type)) +
  geom_bar(stat = 'identity', width = .5, position = 'dodge', color = 'black', alpha = 0.8) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_fill_manual(values = c(brewer.pal(11, 'PuOr')[2], brewer.pal(11, 'PuOr')[9]),
                    labels = c('New repeats', 'Old repeats'), name = 'Repeat type') +
  labs(x="", y="Number of reads per 100 bp") +
  #scale_x_discrete(labels= c('New expansion', 'Old expansion')) +
  ylim(0, 150)


ggsave(plot = read_plot_hyli_hywi,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/read_plot_hyli_hywi.tiff',
       width = 10, height = 10, device = 'tiff')

#------binned Kimura distances-----------------

df_binned = read.delim('/Users/alenakizenko/Documents/PhD/project/pirna/Kimura_binning/LINE_CR1_stats.txt',
                     sep = '\t', header = F)
names(df_binned) = c('Sequence_name', 'Sequence_length', 'Mapped_read-segments', 'Unmapped_read-segments')
df_binned$Normalized_mapped = df_binned$`Mapped_read-segments`/df_binned$Sequence_length

str(df_binned)

read_plot_binned = ggplot(df_binned[-5,], aes(reorder(Sequence_name, Normalized_mapped), Normalized_mapped*100)) +
  geom_bar(stat = 'identity', width = .5, position = 'dodge', color = 'black', alpha = 0.8) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1, angle = 45),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_fill_manual(values = c(brewer.pal(11, 'PuOr')[2], brewer.pal(11, 'PuOr')[9]),
                    labels = c('New repeats', 'Old repeats'), name = 'Repeat type') +
  labs(x="", y="Number of reads per 100 bp") +
  #scale_x_discrete(labels= c('New expansion', 'Old expansion')) +
  ylim(0, 1600)


ggsave(plot = read_plot_hyli_hywi,
       filename = '/Users/alenakizenko/Documents/PhD/project/pictures/read_plot_hyli_hywi.tiff',
       width = 10, height = 10, device = 'tiff')

#read in entire BAM file
bam <- scanBam("/Users/alenakizenko/Documents/PhD/pirna_aep_analysis/hywi/hywi_align_sorted.bam")

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
# expansion
#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#-------------------------------expansion---------------------------------------
#store the mapped positions on the plus and minus strands
exp_neg <- bam_df[bam_df$rname == 'LINE_CR1_expansion' &
                    apply(as.data.frame(bam_df$flag), 1, check_neg),
                  'pos'
]

length(exp)

exp_pos <- bam_df[bam_df$rname == 'LINE_CR1_expansion' &
                    apply(as.data.frame(bam_df$flag), 1, check_pos),
                  'pos'
]
length(rst)

#calculate the densities
exp_neg_density <- density(exp_neg)
exp_pos_density <- density(exp_pos)
#display the negative strand with negative values
exp_neg_density$y <- exp_neg_density$y * -1

exp_final = data.frame('Coordinate_pos' = exp_pos_density$x, 'Coordinate_neg' = exp_neg_density$x,
                       'Density_pos' = exp_pos_density$y, 'Density_neg' = exp_neg_density$y)

exp_pirna_plot = ggplot(exp_final, aes(Coordinate_pos, Density_pos)) +
  geom_line(stat = 'identity', color = brewer.pal(11, 'Spectral')[1]) +
  geom_line(aes(Coordinate_neg, Density_neg), stat = 'identity',
            color = brewer.pal(11, 'Spectral')[10]) +
  geom_hline(yintercept = 0, color = 'gray55') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  #ylim(-8e-07, 8e-07) +
  labs(x="Hydra vulgaris specific expansion", y="Density")

ggsave(plot = exp_pirna_plot,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/hvulg_pirna_exp_hywi.tiff',
       width = 10, height = 10, device = 'tiff')

#----------------------------rest-----------------------------------------------
#store the mapped positions on the plus and minus strands
rest_neg <- bam_df[bam_df$rname == 'LINE_CR1_rest' &
                     apply(as.data.frame(bam_df$flag), 1, check_neg),
                   'pos'
]

rest_pos <- bam_df[bam_df$rname == 'LINE_CR1_rest' &
                     apply(as.data.frame(bam_df$flag), 1, check_pos),
                   'pos'
]

#calculate the densities
rest_neg_density <- density(rest_neg)
rest_pos_density <- density(rest_pos)
#display the negative strand with negative values
rest_neg_density$y <- rest_neg_density$y * -1

rest_final = data.frame('Coordinate_pos' = rest_pos_density$x, 'Coordinate_neg' = rest_neg_density$x,
                        'Density_pos' = rest_pos_density$y, 'Density_neg' = rest_neg_density$y)
rest_pirna_plot = ggplot(rest_final, aes(Coordinate_pos, Density_pos)) +
  geom_line(stat = 'identity', color = brewer.pal(11, 'Spectral')[1]) +
  geom_line(aes(Coordinate_neg, Density_neg), stat = 'identity',
            color = brewer.pal(11, 'Spectral')[10]) +
  geom_hline(yintercept = 0, color = 'gray55') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  # ylim(-8e-07, 8e-07) +
  labs(x="Brown hydra specific expansion", y="Density")

ggsave(plot = rest_pirna_plot,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/hvulg_pirna_rest_hywi.tiff',
       width = 10, height = 10, device = 'tiff')
