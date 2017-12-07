#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

tables <- args[0:(length(args) - 1)]
out.path <- args[length(args)]

IDs <- str_replace(tables, '.+_([A-Z]{3,}\\d{3,}).*', '\\1')
df.list <- map(tables, ~ read_tsv(.x, na = c("-", "NA"))) 
combined <- c()
for (n in seq(length(df.list))){
    df.list[[n]] <- df.list[[n]][1:2] %>% rename('species'=Assemblies)
    df.list[[n]]$sample <- IDs[n]
    df.list[[n]] <- df.list[[n]][c(1,3,2)]
    df.list[[n]][is.na(df.list[[n]])] <- 0
    colnames(df.list[[n]]) <- c('species', 'sample', 'variable')
    print(df.list[[n]])
    combined <- bind_rows(combined, df.list[[n]])
}

combined %>% write_tsv(out.path)




