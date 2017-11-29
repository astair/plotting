#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

tables <- args[0:(length(args) - 1)]
out.path <- args[length(args)]

IDs <- str_replace(tables, '.+/([A-Z]\\w+).\\d+.count.tsv', '\\1')

print(IDs)

df.list <- map(tables, ~ read_tsv(.x, na = c("-", "NA"), col_names = F)) %>% 
  map(function(x) {
    colnames(x)[1] <- 'species'
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
  gather(sample, read_count, -species) %>%
  write_tsv(out.path)