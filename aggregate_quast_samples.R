#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

tables <- args[0:(length(args) - 1)]
out.path <- args[length(args)]

IDs <- str_replace(tables, '.+_([A-Z]{3,}\\d{3,}).*', '\\1')

df.list <- map(tables, ~ read_tsv(.x, na = c("-", "NA"))) %>%
  map(function(x) {
  colnames(x) <- c('species', 'V1') 
  return(x)
  }) %>%
  map( ~ replace_na(.x, list(V1=0)))
  

names(df.list) <- IDs
combined.wide <- df.list[[1]][1]

for (n in seq(length(IDs))){
  col <- df.list[[n]][2]
  colnames(col) <- IDs[n]
  combined.wide <- bind_cols(combined.wide, col)
}

combined.long <- combined.wide %>% 
  gather(sample, V1, -species) %>%
  write_tsv(out.path)




