#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

tables <- args[0:(length(args) - 1)]
out.path <- args[length(args)]

IDs <- str_replace(tables, '.+/([A-Z]\\w+).txt', '\\1')

print(IDs)

df.list <- map(tables, ~ read_tsv(.x, na = c("-", "NA"), col_names = F)) %>% 
  lapply(function(x) {
    colnames(x) <- c('species', 'V2')
    return(x)
  }) 

names(df.list) <- IDs
combined.wide <- df.list[[1]][1]

for (n in seq(length(IDs))){
  col <- df.list[[n]][2]
  colnames(col) <- IDs[n]
  combined.wide <- cbind(combined.wide, col)
}

combined.long <- combined.wide %>%
  gather(sample, V1, -species) %>%
  write_tsv(out.path)