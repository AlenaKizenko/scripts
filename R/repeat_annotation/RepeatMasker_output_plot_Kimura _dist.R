library(ggplot2)
library(reshape2)

#---------------------plotting result file with Kimura substitution levels-------------------

parse_repeat_r1 = read.csv('PATH/TO/INPUT/FILE/landscape.Div.Rname.tab', sep = '\t',
                           stringsAsFactors = F, skip = 1) # reading the file with ending landscape.Div.Rname.tab

colnames(parse_repeat_r1)[4:53] = c(0:49) # naming the columns to bins numbers
parse_repeat_r1$`Repeat family` = paste(parse_repeat_r1$Rclass, parse_repeat_r1$Rfam, sep="/") # merging Class and Family names

# plot all repeat classes

colnames(parse_repeat_r1)[2] = 'Repeat class' # renaming the column

df_repeats_r1_cl = melt(parse_repeat_r1[, -c(1,3,54)], id.vars = 'Repeat class') # reformating the dataframe
df_repeats_r1_cl$perc = c((df_repeats_r1_cl$value*100)/2540000000) # getting percentages - AEP genome size is 254 Mb, NEED TO CHANGE FOR EVERY GENOME

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
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) # plotting

ggsave(plot = plot_r1_plan,
       filename = 'PATH/TO/OUTPUT/plot.tiff',
       width = 10, height = 10, device = 'tiff') # saving the plot