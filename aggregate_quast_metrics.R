#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

tables <- sort(args[0:(length(args) - 1)])
print(tables)
out.path <- args[length(args)]

df.list <- map(tables, ~ read_tsv(.x))
out.table <- df.list[[1]]

for (table in df.list[2:length(df.list)]){
  out.table <- inner_join(out.table, table, by = c('species', 'sample'))
}
variables <- sort(c("#Ns_per_100bp", "#contigs", "#genes", "#indels_per_100kb", "misassemblies", "mismatches_per_100kb", "#operons", "#predicted_genes", "av_contig_read_support", "duplication_ratio", "genome_fraction", "LGA50", "largest_alignment", "largest_contig", "missassembled_contigs_length", "NGA50", "total_aligned_legnth", "total_length_10000bp", "total_length_1000bp", "total_length_50000bp", "total_length"))
names.vec <- c("species", "sample", variables)
colnames(out.table) <- names.vec[1:length(colnames(out.table))]
print(out.table)
write_tsv(out.table, out.path)

