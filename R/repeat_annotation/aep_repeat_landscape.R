library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(Polychrome)
library(wesanderson)

#-----------------------------------1st round-------------------------------------------------------------------------------------------------------

parse_repeat_r1 = read.csv('/Users/alenakizenko/Documents/PhD/aep/aep.final.genome.fa.out.landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F, skip = 1)

names = c(0:49)
colnames(parse_repeat_r1)[4:53] = names
parse_repeat_r1$`Repeat family` = paste(parse_repeat_r1$Rclass, parse_repeat_r1$Rfam, sep="/")

# plot all repeat classes

colnames(parse_repeat_r1)[2] = 'Repeat class'
#parse_repeat_r1 = subset(parse_repeat_r1, parse_repeat_r1$`Repeat class` != 'DNA?')

df_repeats_r1_cl = melt(parse_repeat_r1[, -c(1,3,54)], id.vars = 'Repeat class')
df_repeats_r1_cl$perc = c((df_repeats_r1_cl$value*100)/1250000000) # genome size is 254 Mb

plot_r1_plan = ggplot(data=df_repeats_r1_cl, aes(variable, perc, fill=`Repeat class`))+
  geom_bar(stat="identity", position = position_stack(reverse = FALSE), alpha = 0.8) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
        legend.position = c(0.73, 0.72),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,4,8,9,7,11)]) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4.8))

ggsave(plot = plot_r1_plan,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/plot_r1_aep.tiff',
       width = 10, height = 10, device = 'tiff')

#-----------------------------------2nd round-------------------------------------------------------------------------------------------------------

parse_repeat_r2 = read.csv('/Users/alenakizenko/Documents/PhD/aep/aep.final.genome.fa.masked.out.landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F, skip = 1)

names_r2= c(0:49)
colnames(parse_repeat_r2)[4:53] = names_r2
parse_repeat_r2$`Repeat family` = paste(parse_repeat_r2$Rclass, parse_repeat_r2$Rfam, sep="/")

# plot all repeat classes

colnames(parse_repeat_r2)[2] = 'Repeat class'
parse_repeat_r2 = subset(parse_repeat_r2, parse_repeat_r2$`Repeat class` != 'SINE?')

df_repeats_r2_cl = melt(parse_repeat_r2[, -c(1,3,54)], id.vars = 'Repeat class')
df_repeats_r2_cl$perc = c((df_repeats_r2_cl$value*100)/1250000000)

plot_r2_plan = ggplot(data=df_repeats_r2_cl, aes(variable, perc, fill=`Repeat class`))+
  geom_bar(stat="identity", position = position_stack(reverse = FALSE), alpha = 0.8) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
        legend.position = c(0.73, 0.72),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,4,8,9,7,11)]) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4.8))

ggsave(plot = plot_r2_plan,
       filename = '/Users/alenakizenko/Documents/PhD/pictures/plot_r2_aep.tiff',
       width = 10, height = 10, device = 'tiff')


