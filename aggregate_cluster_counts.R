#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

tables <- args[0:(length(args) - 1)]
out.path <- args[length(args)]

IDs <- str_replace(tables, '.+/([A-Z]\\w+).\\d+.count.tsv', '\\1')
clusters <- str_replace(tables, '.+/[A-Z]\\w+.(\\d+).count.tsv', '\\1')

print(IDs)
print(clusters)

df.list <- map(tables, ~ read_tsv(.x, na = c("-", "NA"), col_names = F)) 

combined.df <- c()
for (n in seq(length(df.list))){
    x <- df.list[[n]]
    if (length(colnames(x)) > 0){
        x$sample <- IDs[n]
        x$cluster <- clusters[n]
        colnames(x) <- c('species', 'read_count', 'read_abundance', 'sample', 'cluster')
        combined.df <- bind_rows(combined.df, x)
    }
}

combined.df %>% write_tsv(out.path)


