library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(Polychrome)
library(wesanderson)

#-------------------------------------hidra viridissima--------------------------------------------------------------------------------------------
#-----------------------------------1st round-------------------------------------------------------------------------------------------------------

parse_repeat_r1 = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/round_1/parse_rep_land_rma1/hvir_genome_hm2_250116_renamed.fa.align.landscape.Div.Rname.tab', sep = '\t',
                            stringsAsFactors = F)
    
names = c(0:49)
colnames(parse_repeat_r1)[4:53] = names
parse_repeat_r1$`Repeat family` = paste(parse_repeat_r1$Rclass, parse_repeat_r1$Rfam, sep="/")

    
df_repeats_r1_fam = melt(parse_repeat_r1[, -c(1,2,3)], id.vars = 'Repeat family')
df_repeats_r1_fam$perc = c((df_repeats_r1_fam$value*100)/254000000) # genome size is 254 Mb

df_repeats_r1_fam = subset(df_repeats_r1_fam, df_repeats_r1_fam$value > 0)

# plot all repeat families
    
ggplot(data=df_repeats_r1_fam, aes(variable, perc, fill=`Repeat family`))+
geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0),
        legend.position = c(0.73, 0.72),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
scale_fill_manual(values = c("DNA/hAT-hATm" = "#FFDB6D", "DNA/TcMar-Tc1" = "#C4961A", "Unknown/Unknown" = "#999999", "DNA/Crypton-H" = "#293352",
                            "DNA/P" = "#D16103", "DNA/Maverick" = "#C3D7A4", "DNA/CMC-Chapaev-3" = "#52854C", "DNA/CMC-Transib" = "#4E84C4",
                            "LINE/CR1" = "#118E53", "LINE/Penelope" = "#E69F00", "DNA/CMC-Chapaev" = "#56B4E9", "DNA/TcMar-Fot1" = "#009E73",
                            "DNA/CMC-EnSpm" = "#F0E442", "DNA/Sola-2" = "#0072B2", "LINE/L2" = "#D55E00", "DNA/TcMar-Mariner" = "#CC79A7",
                            "DNA/hAT-Ac" = "#F1AB5D", "DNA/hAT-hATx" = "#C79F1A", "DNA/Zator" = "#14A1C0", "DNA/PIF-Harbinger" = "#203312",
                            "DNA/hAT-Tip100" =  "#A19005", "DNA/TcMar-Pogo" = "#C3F7A2", "DNA/hAT-hAT5" = "#3285DC", "DNA/Academ-1" = "#1E44F4",
                            "DNA/Ginger-1" = "#666666", "DNA/hAT-Tag1" = "#D19A52", "RC/Helitron" = "#52B9B9", "DNA/PiggyBac" = "#229A13",
                            "DNA/MULE-MuDR" = "#DDE147", "LTR/Gypsy" = "#5074D2", "DNA/Sola-1" = "#D15A61", "LINE/L1-Tx1" = "#AD7917",
                            "LTR/Pao" = "#D0E541", "DNA/Merlin" = "#1174A2", "DNA/hAT-hATw" = "#D25D01", "DNA/TcMar-Stowaway" = "#AC76A6",
                            "DNA/Kolobok-Hydra" = "#A1DD6F", "LINE/R2" = "#A4465A", "DNA/Ginger" = "#F2EADA", "DNA/hAT-Charlie" = "#143454",
                            "DNA/hAT-hobo" = "#CCCCCC", "DNA/TcMar-Tc2" = "#A77D10", "DNA/DNA" = "#46A4E7", "LINE/L1" = "#F4EDCA",
                            "DNA/hAT-Blackjack" = "#D4ADAC", "DNA/PIF-ISL2EU" = "#1D2543", "Other/DNA_virus" = "#E57E11", "SINE/tRNA" = "#776F72",
                            "DNA/TcMar-ISRm11" = "#C2562A", "LINE/RTE-X" = "#338D57", "DNA/Sola-3" = "#37A4D8", "DNA/hAT" = "#F9DF6A",
                            "DNA/hAT-Pegasus" = "#EFAC2A", "DNA/hAT-hAT19" = "#C5962C", "DNA/Academ-2" = "#015EA3"
))


# plot only CR1 repeats

df_repeats_r1_fam_cr1 = subset(df_repeats_r1_fam, `Repeat family` == "LINE/CR1")

ggplot(data=df_repeats_r1_fam_cr1, aes(variable, perc))+
  geom_bar(stat="identity", position = position_stack(reverse = FALSE), fill="#118E53") +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0),
        legend.position = "none",
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50))


# plot all repeat classes

colnames(parse_repeat_r1)[2] = 'Repeat class'

df_repeats_r1_cl = melt(parse_repeat_r1[, -c(1,3,54)], id.vars = 'Repeat class')
df_repeats_r1_cl$perc = c((df_repeats_r1_cl$value*100)/254000000) # genome size is 254 Mb

plot_r1_hvir = ggplot(data=df_repeats_r1_cl, aes(variable, perc, fill=`Repeat class`))+
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

ggsave(plot = plot_r1_hvir,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/plot_r1_hvir.tiff',
       width = 10, height = 10, device = 'tiff')

#-------------------------------pie plot round 1------------------------------------------------------------------------------------

df_pie_r1_fam = df_repeats_r1_fam[, c(1,3)]

df_pie_r1_fam = aggregate(df_pie_r1_fam$value, by = list(df_pie_r1_fam$`Repeat family`), FUN = sum)
colnames(df_pie_r1_fam) = c('Repeat family', 'Size')

df_pie_r1_fam$perc = c((df_pie_r1_fam$Size*100)/254000000)

nonrep_nt_r1 = 254000000 - sum(df_pie_r1_fam$Size)

nonrep_pers_r1 = nonrep_nt_r1*100/254000000

df_pie_r1_fam = rbind(df_pie_r1_fam, c('Non-repetitive', nonrep_nt_r1, nonrep_pers_r1))

df_pie_r1_fam$Size = as.numeric(df_pie_r1_fam$Size)
df_pie_r1_fam$perc = as.numeric(df_pie_r1_fam$perc)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

df_pie_r1_fam$perc = round(df_pie_r1_fam$perc, 1)

df_pie_r1_fam <- df_pie_r1_fam %>%
  arrange(desc(`Repeat family`)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)


ggplot(df_pie_r1_fam, aes(x='', y=perc, fill = `Repeat family`)) +
  geom_bar(stat = 'identity', colour = 'black') +
  coord_polar('y', start = 0) +
  scale_fill_manual(values = c("DNA/hAT-hATm" = "#FFDB6D", "DNA/TcMar-Tc1" = "#C4961A", "Unknown/Unknown" = "#999999", "DNA/Crypton-H" = "#293352",
                               "DNA/P" = "#D16103", "DNA/Maverick" = "#C3D7A4", "DNA/CMC-Chapaev-3" = "#52854C", "DNA/CMC-Transib" = "#4E84C4",
                               "LINE/CR1" = "#118E53", "LINE/Penelope" = "#E69F00", "DNA/CMC-Chapaev" = "#56B4E9", "DNA/TcMar-Fot1" = "#009E73",
                               "DNA/CMC-EnSpm" = "#F0E442", "DNA/Sola-2" = "#0072B2", "LINE/L2" = "#D55E00", "DNA/TcMar-Mariner" = "#CC79A7",
                               "DNA/hAT-Ac" = "#F1AB5D", "DNA/hAT-hATx" = "#C79F1A", "DNA/Zator" = "#14A1C0", "DNA/PIF-Harbinger" = "#203312",
                               "DNA/hAT-Tip100" =  "#A19005", "DNA/TcMar-Pogo" = "#C3F7A2", "DNA/hAT-hAT5" = "#3285DC", "DNA/Academ-1" = "#1E44F4",
                               "DNA/Ginger-1" = "#666666", "DNA/hAT-Tag1" = "#D19A52", "RC/Helitron" = "#52B9B9", "DNA/PiggyBac" = "#229A13",
                               "DNA/MULE-MuDR" = "#DDE147", "LTR/Gypsy" = "#5074D2", "DNA/Sola-1" = "#D15A61", "LINE/L1-Tx1" = "#AD7917",
                               "LTR/Pao" = "#D0E541", "DNA/Merlin" = "#1174A2", "DNA/hAT-hATw" = "#D25D01", "DNA/TcMar-Stowaway" = "#AC76A6",
                               "DNA/Kolobok-Hydra" = "#A1DD6F", "LINE/R2" = "#A4465A", "DNA/Ginger" = "#F2EADA", "DNA/hAT-Charlie" = "#143454",
                               "DNA/hAT-hobo" = "#CCCCCC", "DNA/TcMar-Tc2" = "#A77D10", "DNA/DNA" = "#46A4E7", "LINE/L1" = "#F4EDCA",
                               "DNA/hAT-Blackjack" = "#D4ADAC", "DNA/PIF-ISL2EU" = "#1D2543", "Other/DNA_virus" = "#E57E11", "SINE/tRNA" = "#776F72",
                               "DNA/TcMar-ISRm11" = "#C2562A", "LINE/RTE-X" = "#338D57", "DNA/Sola-3" = "#37A4D8", "DNA/hAT" = "#F9DF6A",
                               "DNA/hAT-Pegasus" = "#EFAC2A", "DNA/hAT-hAT19" = "#C5962C", "DNA/Academ-2" = "#015EA3", "Non-repetitive" = "#DDDDDD"
  )) +
  blank_theme +
  theme(axis.text.x=element_blank(),
        legend.position = 'none') +
  geom_text(aes(y = lab.ypos, label = ifelse(perc > 1, perc, "")), color = "black", size = 6) +
  geom_text(aes(label = ifelse(perc > 10, `Repeat family`, "")), color = "black", size = 6)
  

#-----------------------------------2nd round-------------------------------------------------------------------------------------------------------

parse_repeat_r2 = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/round_2_after_hardmasking/parse_rep_land_rma2/hvir_genome_hm2_250116_renamed.fa.masked.align.landscape.Div.Rname.tab', sep = '\t',
                        stringsAsFactors = F)

names_r2= c(0:49)
colnames(parse_repeat_r2)[4:53] = names_r2
parse_repeat_r2$`Repeat family` = paste(parse_repeat_r2$Rclass, parse_repeat_r2$Rfam, sep="/")


df_repeats_r2 = melt(parse_repeat_r2[, -c(1,2,3)], id.vars = 'Repeat family')
df_repeats_r2$perc = c((df_repeats_r2$value*100)/254000000) # genome size is 254 Mb
df_repeats_r2 = subset(df_repeats_r2, df_repeats_r2$value > 0)

ggplot(data=df_repeats_r2, aes(variable, perc, fill=`Repeat family`))+
  geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
        legend.position = c(0.73, 0.72),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  scale_fill_manual(values = c("DNA/TcMar-Mariner" = "#CC79A7", "Unknown/Unknown" = "#999999", "DNA/Sola-2" = "#0072B2", "LINE/L1-Tx1" = "#AD7917",
                               "DNA/CMC-Chapaev" = "#56B4E9", "DNA/Zator" = "#14A1C0", "DNA/Maverick" = "#C3D7A4", "DNA/hAT-hATm" = "#FFDB6D",
                               "DNA/CMC-Chapaev-3" = "#52854C", "DNA/TcMar-Stowaway" = "#AC76A6", "LTR/Gypsy" = "#5074D2", "DNA/CMC-Transib" = "#4E84C4",
                               "DNA/CMC-EnSpm" = "#F0E442", "LINE/L2" = "#D55E00", "DNA/TcMar-Tc1" = "#C4961A", "DNA/TcMar-Tc2" = "#A77D10",
                               "DNA/hAT-hATx" = "#C79F1A", "DNA/Kolobok-Hydra" = "#A1DD6F", "Other/DNA_virus" = "#E57E11", "DNA/PIF-ISL2EU" = "#1D2543",
                               "DNA/hAT-Tip100" =  "#A19005", "DNA/P" = "#D16103", "DNA/DNA" = "#46A4E7", "DNA/Sola-1" = "#D15A61",
                               "DNA/hAT-Ac" = "#F1AB5D", "DNA/TcMar-Mariner" = "#CC79A7", "DNA/TcMar-Pogo" = "#C3F7A2", "DNA/PIF-Harbinger" = "#203312",
                               "LINE/R2" = "#A4465A", "DNA/PiggyBac" = "#229A13", "DNA/hAT-Tag1" = "#D19A52", "LINE/CR1" = "#118E53", "LTR/Pao" = "#D0E541",
                               "DNA/hAT" = "#F9DF6A", "RC/Helitron" = "#52B9B9", "DNA/TcMar-Fot1" = "#009E73"))

# plot CR1 family

df_repeats_r2_cr1 = subset(df_repeats_r2, `Repeat family` == "LINE/CR1")

ggplot(data=df_repeats_r2_cr1, aes(variable, perc))+
  geom_bar(stat="identity", position = position_stack(reverse = FALSE), fill="#118E53") +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion (%)")+
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0),
        legend.position = "none",
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50))


# plot all repeat classes

colnames(parse_repeat_r2)[2] = 'Repeat class'

df_repeats_r2_cl = melt(parse_repeat_r2[, -c(1,3,54)], id.vars = 'Repeat class')
df_repeats_r2_cl$perc = c((df_repeats_r2_cl$value*100)/254000000) # genome size is 254 Mb

plot_r2_hvir = ggplot(data=df_repeats_r2_cl, aes(variable, perc, fill=`Repeat class`))+
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
  scale_fill_manual(name = 'Repeat class', values = brewer.pal(11, 'Spectral')[c(1,2,4,8,9,11)]) +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4.8))

ggsave(plot = plot_r2_hvir,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/plot_r2_hvir.tiff',
       width = 10, height = 10, device = 'tiff')

#----------------------------------------pie plot round 2--------------------------------------------------------------------------------------------------

df_pie_r2 = df_repeats_r2[, c(1,3)]

df_pie_r2 = aggregate(df_pie_r2$value, by = list(df_pie_r2$`Repeat family`), FUN = sum)
colnames(df_pie_r2) = c('Repeat family', 'Size')

df_pie_r2$perc = c((df_pie_r2$Size*100)/254000000)
sum(df_pie_r2$Size)

df_pie_r2 = rbind(df_pie_r2, c('Non-repetitive', '247032194', '97.25677'))

df_pie_r2$Size = as.numeric(df_pie_r2$Size)
df_pie_r2$perc = as.numeric(df_pie_r2$perc)


df_pie_r2$perc = round(df_pie_r2$perc, 1)

df_pie_r2 <- df_pie_r2 %>%
  arrange(desc(`Repeat family`)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)


ggplot(df_pie_r2, aes(x='', y=perc, fill = `Repeat family`)) +
  geom_bar(stat = 'identity', colour = 'black') +
  coord_polar('y', start = 0) +
  scale_fill_manual(values = c("DNA/TcMar-Mariner" = "#CC79A7", "Unknown/Unknown" = "#999999", "DNA/Sola-2" = "#0072B2", "LINE/L1-Tx1" = "#AD7917",
                               "DNA/CMC-Chapaev" = "#56B4E9", "DNA/Zator" = "#14A1C0", "DNA/Maverick" = "#C3D7A4", "DNA/hAT-hATm" = "#FFDB6D",
                               "DNA/CMC-Chapaev-3" = "#52854C", "DNA/TcMar-Stowaway" = "#AC76A6", "LTR/Gypsy" = "#5074D2", "DNA/CMC-Transib" = "#4E84C4",
                               "DNA/CMC-EnSpm" = "#F0E442", "LINE/L2" = "#D55E00", "DNA/TcMar-Tc1" = "#C4961A", "DNA/TcMar-Tc2" = "#A77D10",
                               "DNA/hAT-hATx" = "#C79F1A", "DNA/Kolobok-Hydra" = "#A1DD6F", "Other/DNA_virus" = "#E57E11", "DNA/PIF-ISL2EU" = "#1D2543",
                               "DNA/hAT-Tip100" =  "#A19005", "DNA/P" = "#D16103", "DNA/DNA" = "#46A4E7", "DNA/Sola-1" = "#D15A61",
                               "DNA/hAT-Ac" = "#F1AB5D", "DNA/TcMar-Mariner" = "#CC79A7", "DNA/TcMar-Pogo" = "#C3F7A2", "DNA/PIF-Harbinger" = "#203312",
                               "LINE/R2" = "#A4465A", "DNA/PiggyBac" = "#229A13", "DNA/hAT-Tag1" = "#D19A52", "LINE/CR1" = "#118E53", "LTR/Pao" = "#D0E541",
                               "DNA/hAT" = "#F9DF6A", "RC/Helitron" = "#52B9B9", "DNA/TcMar-Fot1" = "#009E73", "Non-repetitive" = "#DDDDDD"
  )) +
  blank_theme +
  theme(axis.text.x=element_blank(),
        legend.position = 'none') +
  geom_text(aes(y = lab.ypos, label = ifelse(perc > 1, perc, "")), color = "black", size = 6) +
  geom_text(aes(y = 9, label = ifelse(perc > 10, `Repeat family`, "")), color = "black", size = 6)


#--------------------repeatcraft--------------------------------------------------------------------------------------------

rptcrft_df = read.csv('/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeatcraft/green_masked_rptcrf_out.txt', sep = '\t', stringsAsFactors = F)
str(rptcrft_df)
colnames(rptcrft_df) = c('Repeat class', 'Before merge', 'After merge')
rptcrft_df = melt(rptcrft_df, id.vars = 'Repeat class')

ggplot(rptcrft_df, aes(x = `Repeat class`, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, vjust = 1, hjust=1, angle = 60),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,3)],
                    name="Number of repeats",
                    # breaks=c("TP53 mutant", "TP53 wild type"),
                    labels=c("Before merge", "After merge")) +
  ylab('Number of repeats') +
  scale_y_continuous(breaks=seq(0, 15000, by = 2000))


#---------------------------------------------------plot distribution of repeats in first and second round------------------------------------------------------------

r1_agg = aggregate(x = df_repeats_r1_cl$value, by = list(df_repeats_r1_cl$`Repeat class`), FUN = sum)
colnames(r1_agg) = c('Repeat class', 'DNA amount')
r1_agg$Round = c(rep('Standard masking', each =7))
r2_agg = aggregate(x = df_repeats_r2_cl$value, by = list(df_repeats_r2_cl$`Repeat class`), FUN = sum)
colnames(r2_agg) = c('Repeat class', 'DNA amount')
r2_agg = rbind(r2_agg, c('SINE', '0'))
r2_agg$Round = c(rep('Double masking', each = 7))         
r12_agg = rbind(r1_agg, r2_agg)
str(r12_agg)
r12_agg$`DNA amount` = as.numeric(r12_agg$`DNA amount`)
r12_agg$Round = as.factor(r12_agg$Round)
r12_agg = r12_agg[order(r12_agg$`DNA amount`),]
r12_agg$perc = (r12_agg$`DNA amount`*100)/254000000
positions = c('DNA', 'Unknown', 'LTR', 'LINE', 'RC', 'Other', 'SINE')
r12_agg$Round <- ordered(r12_agg$Round, levels = c('Standard masking', 'Double masking'))

dist_r12 = ggplot(r12_agg, aes(`Repeat class`, perc, fill = Round)) +
  geom_bar(stat = 'identity', width = .5, alpha= .8) +
  scale_x_discrete(limits = positions) +
  facet_grid(cols = vars(Round)) +
  theme_bw() +
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(3,2)]) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 18, colour = "black", angle = 0),
        legend.position = 'none',
        axis.title.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 1, hjust=1, angle = 60, color = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title=element_text(size=18)) +
  ylab("Genome Proportion (%)") +
  ggtitle('Distribution of repeat families after standard and double masking') +
  scale_y_continuous(breaks=seq(0, 30, 5), limits=c(0, 30))

ggsave(plot = dist_r12,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeatcraft/dist_r12_hvir.jpg',
       width = 10, height = 10, device = 'jpeg')

ggsave(plot = dist_r12,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeatcraft/dist_r12_hvir.tiff',
       width = 10, height = 10, device = 'tiff')

#---------combined 1st and 2nd rounds repeat landscape------------------------

title <- ggdraw() + 
  draw_label(
    "Repeat landscape of Hydra viridissima genome",
    fontface = 'bold', size = 20,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(-10, 0, 0, 370)
  )

plot_row <- plot_grid(plot_r1_hvir, plot_r2_hvir, labels = c('                        Standard masking',
                                                             '                        Double masking'),
                      label_size = 20,
                      #label_x = 0, label_y = 0,
                      hjust = -0.28, vjust = -0.5,
                      rel_widths = c(1, 1), rel_heights = c(1,1))

combined_plots = plot_grid(title, plot_row,
                           ncol = 1,
                           # rel_heights values control vertical title margins
                           rel_heights = c(0.1, 1))

ggsave(plot = combined_plots,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeatcraft/combined_plots_repeatlandscape_hvir.jpeg',
       width = 15, height = 10, device = 'jpeg')

ggsave(plot = combined_plots,
       filename = '/Users/alenakizenko/Documents/PhD/green_hydra_repeat_ann/repeatcraft/combined_plots_repeatlandscape_hvir.tiff',
       width = 15, height = 10, device = 'tiff')

