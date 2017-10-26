#!/usr/bin/env Rscript

library(tidyverse)
library(viridis)

args = commandArgs(trailingOnly=TRUE)

RA <- read.delim(args[1])
GF <- read.delim(args[2])
out.file <- args[3]

funky <- c('#B71C1C', '#0D47A1', '#388E3C', '#F9A825', '#E64A19', '#1976D2')
black <- rep('black', 50)

GF.sorted <- GF %>% 
  group_by(sample) %>% 
  arrange(V1, .by_group=T) 

comb <- inner_join(GF.sorted, RA, by=c('species', 'sample'))

p = ggplot(comb, aes(x=log10(V1.y), y=V1.x, color=species)) + 
  geom_point(size=2, alpha=0.8) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  # scale_color_manual(values=c(funky, black)) + 
  scale_color_viridis(discrete=T, option='magma') + 
  theme(legend.position="none") + 
  labs(x='log(Relative Abundance)', y='Genome Faction [%]')

pdf(out.file, width=6, height=4)
p
dev.off()