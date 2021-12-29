library(ggplot2)
library(circlize)
library(RColorBrewer)
library(viridis)
library(plyr)
library(reshape2)
library(tidyr)
library(ggsignif)
library(gg.gap)
library(cowplot)
library(tidyverse)
library(cowplot)
library(pafr)
library(dplyr)


chr_dov = read.delim2('/Users/alenakizenko/Documents/PhD/project/hydra_vulgaris_chrom/hvul_chr_dovtail_align.paf',
                      sep = '\t', header = F)

