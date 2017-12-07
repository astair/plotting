library(tidyverse)

setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')

nnames <- c('species', 'seq_num', 'assembly_len', 'mean_len', 'longest_contig', 'shortest_contig', 'GC_cont', 'N_cont', 'N50', 'L50')

genome_stats_complex <- read_tsv(
    'nile/data/specs/complex_genomes_stats.txt', col_names = F) %>% 
    rename_all( ~nnames) 
genome_stats_low <- read_tsv(
    'nile/data/specs/low_genomes_stats.txt', col_names = F) %>%
    rename_all( ~nnames) 
genome_stats_med <- read_tsv(
    'nile/data/specs/med_genomes_stats.txt', col_names = F) %>%
    rename_all( ~nnames) 


### GENOME STATS AND REALTIVE ABUNDANCES ###

seq_counts_complex <- read_tsv('results/counts_complex.tsv') 
seq_counts_med <- read_tsv('results/counts_med.tsv') 
seq_counts_low <- read_tsv('results/counts_low.tsv') 

relabu_complex <- read_tsv('results/abundance_complex.tsv') %>% 
    rename('rel_abundance'=V1)
relabu_med <- read_tsv('results/abundance_med.tsv') %>% 
    rename('rel_abundance'=V1)
relabu_low <- read_tsv('results/abundance_low.tsv') %>% 
    rename('rel_abundance'=V1)

info_complex <- inner_join(genome_stats_complex, seq_counts_complex) %>% 
    inner_join(relabu_complex) %>% 
    mutate(genome_coverage = (read_count * 77) / assembly_len) 
    # %>% write_tsv('results/info_complex.tsv') 
info_med <- inner_join(genome_stats_med, seq_counts_med) %>% 
    inner_join(relabu_med) %>% 
    mutate(genome_coverage = (read_count * 77) / assembly_len) 
    # %>% write_tsv('results/info_med.tsv') 
info_low <- inner_join(genome_stats_low, seq_counts_low) %>% 
    inner_join(relabu_low) %>% 
    mutate(genome_coverage = (read_count * 77) / assembly_len) 
    # %>% write_tsv('results/info_low.tsv') 

info_all <- bind_rows(complex=info_complex, medium=info_med, low=info_low, .id = 'complexity')


### SINGLE SAMPLE ASSEMBLY ###

megahit_500b_med <- read_tsv(
    'results/final/megahit_medium_500b.tsv') 
megahit_10kb_med <- read_tsv(
    'results/final/megahit_medium_10kb.tsv') 
megahit_500b_complex <- read_tsv(
    'results/final/megahit_complex_500b.tsv') 
megahit_10kb_complex <- read_tsv(
    'results/final/megahit_complex_10kb.tsv') 
megahit_500b_low <- read_tsv(
    'results/final/megahit_low_500b.tsv') 
megahit_10kb_low <- read_tsv(
    'results/final/megahit_low_10kb.tsv') 

spades_500b_med <- read_tsv(
    'results/final/spades_medium_500b.tsv') 
spades_10kb_med <- read_tsv(
    'results/final/spades_medium_10kb.tsv') 
spades_500b_complex <- read_tsv(
    'results/final/spades_complex_500b.tsv') 
spades_10kb_complex <- read_tsv(
    'results/final/spades_complex_10kb.tsv') 
spades_500b_low <- read_tsv(
    'results/final/spades_low_500b.tsv') 
spades_10kb_low <- read_tsv(
    'results/final/spades_low_10kb.tsv') 

med_500b <- bind_rows(
    MEGAHIT=megahit_500b_med, SPADES=spades_500b_med, .id = 'assembler')
med_10kb <- bind_rows(
    MEGAHIT=megahit_10kb_med, SPADES=spades_10kb_med, .id = 'assembler')

complex_500b <- bind_rows(
    MEGAHIT=megahit_500b_complex, SPADES=spades_500b_complex, .id = 'assembler')
complex_10kb <- bind_rows(
    MEGAHIT=megahit_10kb_complex, SPADES=spades_10kb_complex, .id = 'assembler')

low_500b <- bind_rows(
    MEGAHIT=megahit_500b_low, SPADES=spades_500b_low, .id = 'assembler')
low_10kb <- bind_rows(
    MEGAHIT=megahit_10kb_low, SPADES=spades_10kb_low, .id = 'assembler')

med <- bind_rows('500b'=med_500b, '10kb'=med_10kb, .id = 'cutoff')
complex <- bind_rows('500b'=complex_500b, '10kb'=complex_10kb, .id = 'cutoff')
low <- bind_rows('500b'=low_500b, '10kb'=low_10kb, .id = 'cutoff')

ss_assembly <- bind_rows(medium=med, complex=complex, low=low, .id = 'complexity')
ss_assembly_info <- inner_join(ss_assembly, info_all)



### BINNERS ###

concoct_500b_med <- read_tsv(
    'results/final/concoct_medium_500b.tsv') 
concoct_10kb_med <- read_tsv(
    'results/final/concoct_medium_10kb.tsv') 
lsa_500b_med <- read_tsv('results/final/lsa_medium_55_500b.tsv') 
lsa_10kb_med <- read_tsv('results/final/lsa_medium_55_10kb.tsv')
lsa_500b_33_med <- read_tsv('results/final/lsa_medium_33_500b.tsv') 
lsa_10kb_33_med <- read_tsv('results/final/lsa_medium_33_10kb.tsv')
# lsa_500b_15_med <- read_tsv('results/final/lsa_medium_15_500b.tsv') 
# lsa_10kb_15_med <- read_tsv('results/final/lsa_medium_15_10kb.tsv')

concoct <- bind_rows('500b'=concoct_500b_med, '10kb'=concoct_10kb_med, .id = 'cutoff')
lsa <- bind_rows('500b'=lsa_500b_med, '10kb'=lsa_10kb_med, .id = 'cutoff')
lsa_33 <- bind_rows('500b'=lsa_500b_33_med, '10kb'=lsa_10kb_33_med, .id = 'cutoff')
# lsa_15 <- bind_rows('500b'=lsa_500b_15_med, '10kb'=lsa_10kb_15_med, .id = 'cutoff')

binned_assembly <- bind_rows('LSA_55'=lsa, 'LSA_33'=lsa_33, 'CONCOCT'=concoct, .id = 'assembler')
binned_assembly$complexity <- 'medium'

binned_assembly_info <- inner_join(binned_assembly, info_all, .by=c(sample, species))

the_total_tbl <- bind_rows(ss_assembly_info, binned_assembly_info) %>% rename('total_aligned_length'=total_aligned_legnth)


### Trying to calculate putity or something ###
concoct_counts_medium <- read_tsv('results/final/concoct_medium_10_count.tsv')
lsa_counts_medium <- read_tsv('results/final/lsa_medium_55_count.tsv')

binners_counts <- bind_rows(CONCOCT=concoct_counts_medium, LSA=lsa_counts_medium, .id='assembler') %>% rename('read_count_cluster'=read_count, 'read_abundance_cluster'=read_abundance)
binners_counts$complexity <- 'medium'
binners_count_info <- inner_join(binners_counts, info_all) %>% group_by(cluster, sample, assembler) %>% 
  mutate(
    purity=max(read_abundance_cluster), 
    cluster_owner=c(species[read_abundance_cluster==purity])[1], 
    rel_abundance_cluster_owner=rel_abundance[species==cluster_owner]
    ) %>% 
  ungroup() %>%
  group_by(species, sample, assembler) %>%
  mutate(
    max_rel_abundance=max(read_abundance_cluster)
  )

### FIDDELING WITH STRAINS ### 
strains_tbl <- read_tsv('nile/data/specs/strains_mapping_2.txt')
ani_tbl <- read_tsv('results/perc_ids.tab', skip=2) %>% gather(key=representative, value=ANI, -X1) %>% rename('species'=X1) %>% print() 