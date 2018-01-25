library(tidyverse)
library(stringr)

setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')

nnames <- c('species', 'seq_num', 'assembly_len', 'mean_len', 'longest_contig', 'shortest_contig', 'GC_cont', 'N_cont', 'N50', 'L50')

genome_stats_complex <- read_tsv(
    'nile/data/sim/complex_genomes_stats.txt', col_names = F) %>% 
    rename_all( ~nnames) 
genome_stats_low <- read_tsv(
    'nile/data/sim/low_genomes_stats.txt', col_names = F) %>%
    rename_all( ~nnames) 
genome_stats_med <- read_tsv(
    'nile/data/sim/med_genomes_stats.txt', col_names = F) %>%
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
    'results/new/megahit_medium_500b.tsv') 
megahit_10kb_med <- read_tsv(
    'results/new/megahit_medium_10kb.tsv') 
megahit_500b_complex <- read_tsv(
    'results/new/megahit_complex_500b.tsv') 
megahit_10kb_complex <- read_tsv(
    'results/new/megahit_complex_10kb.tsv') 
megahit_500b_low <- read_tsv(
    'results/new/megahit_low_500b.tsv') 
megahit_10kb_low <- read_tsv(
    'results/new/megahit_low_10kb.tsv') 

spades_500b_med <- read_tsv(
    'results/new/spades_medium_500b.tsv') 
spades_10kb_med <- read_tsv(
    'results/new/spades_medium_10kb.tsv') 
spades_500b_complex <- read_tsv(
    'results/new/spades_complex_500b.tsv') 
spades_10kb_complex <- read_tsv(
    'results/new/spades_complex_10kb.tsv') 
spades_500b_low <- read_tsv(
    'results/new/spades_low_500b.tsv') 
spades_10kb_low <- read_tsv(
    'results/new/spades_low_10kb.tsv') 

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

ss_assembly <- bind_rows(medium=med, complex=complex, low=low, .id = 'complexity') %>% 
    rename('misassembled_contig_length'=misassemblies, 'misassemblies'=mismatches_per_100kb, 'mismatches_per_100kb'=missassembled_contigs_length)

ss_assembly_info <- inner_join(ss_assembly, info_all)



### BINNERS ###

concoct_500b_med <- read_tsv(
    'results/new/concoct_medium_500b.tsv') 
concoct_10kb_med <- read_tsv(
    'results/new/concoct_medium_10kb.tsv') 
lsa_500b_med <- read_tsv('results/new/lsa_medium_55_500b.tsv') 
lsa_10kb_med <- read_tsv('results/new/lsa_medium_55_10kb.tsv')

lsa_500b_33_med <- read_tsv('results/final/lsa_medium_33_500b.tsv') %>% 
    rename('misassembled_contig_length'=misassemblies, 'misassemblies'=mismatches_per_100kb, 'mismatches_per_100kb'=missassembled_contigs_length)
lsa_10kb_33_med <- read_tsv('results/final/lsa_medium_33_10kb.tsv') %>% 
    rename('misassembled_contig_length'=misassemblies, 'misassemblies'=mismatches_per_100kb, 'mismatches_per_100kb'=missassembled_contigs_length)

lsa_500b_15_med <- read_tsv('results/new/lsa_medium_15_500b.tsv') %>% 
    rename('misassembled_contig_length'=misassemblies, 'misassemblies'=mismatches_per_100kb, 'mismatches_per_100kb'=missassembled_contigs_length)
lsa_10kb_15_med <- read_tsv('results/new/lsa_medium_15_10kb.tsv') %>% 
    rename('misassembled_contig_length'=misassemblies, 'misassemblies'=mismatches_per_100kb, 'mismatches_per_100kb'=missassembled_contigs_length)

concoct <- bind_rows('500b'=concoct_500b_med, '10kb'=concoct_10kb_med, .id = 'cutoff')
lsa <- bind_rows('500b'=lsa_500b_med, '10kb'=lsa_10kb_med, .id = 'cutoff')
lsa_33 <- bind_rows('500b'=lsa_500b_33_med, '10kb'=lsa_10kb_33_med, .id = 'cutoff')
lsa_15 <- bind_rows('500b'=lsa_500b_15_med, '10kb'=lsa_10kb_15_med, .id = 'cutoff')

binned_assembly <- bind_rows('LSA_55'=lsa, 'LSA_33'=lsa_33, 'LSA_15'=lsa_15, 'CONCOCT'=concoct, .id = 'assembler')
binned_assembly$complexity <- 'medium'


unique(binned_assembly %>% filter(assembler == 'LSA_15') %>% select(sample))
binned_assembly %>% filter(assembler == 'LSA_15')

unique(binned_assembly %>% filter(assembler == 'LSA_33') %>% select(sample))
binned_assembly %>% filter(assembler == 'LSA_33')

binned_assembly_info <- inner_join(binned_assembly, info_all, .by=c(sample, species))

the_total_tbl <- bind_rows(ss_assembly_info, binned_assembly_info) %>% rename('total_aligned_length'=total_aligned_legnth)


### Trying to calculate putity or something ###
concoct_counts_medium <- read_tsv('results/new/concoct_medium_10_count.tsv')
lsa_counts_medium <- read_tsv('results/new/lsa_medium_55_count.tsv')
lsa_33_counts_medium <- read_tsv('results/new/lsa_medium_33_count.tsv')
lsa_15_counts_medium <- read_tsv('results/new/lsa_medium_15_count.tsv')


binners_counts <- bind_rows(CONCOCT=concoct_counts_medium, LSA_55=lsa_counts_medium, LSA_33=lsa_33_counts_medium, LSA_15=lsa_15_counts_medium, .id='assembler') %>% rename('read_count_cluster'=read_count, 'read_abundance_cluster'=read_abundance)
binners_counts$complexity <- 'medium'
binners_count_info <- inner_join(binners_counts, info_all) %>% 
    group_by(cluster, sample, assembler) %>% 
    mutate(
        purity=max(read_abundance_cluster), 
        cluster_owner=c(ifelse(
            read_abundance_cluster > 0.5, species[read_abundance_cluster==purity][1], NA)),
        rel_abundance_cluster_owner=rel_abundance[species==cluster_owner]
    ) %>% 
    ungroup() %>%
    group_by(species, sample, assembler) %>%
    mutate(max_rel_abundance=max(read_abundance_cluster)
  )

binners_count_across <- binners_count_info %>% 
    select(species, cluster, read_count_cluster, assembler) %>% 
    group_by(cluster, species, assembler) %>% 
    summarise(count=sum(read_count_cluster)) %>% 
    ungroup() %>% 
    group_by(cluster, assembler) %>% 
    mutate(abundance=count/sum(as.numeric(count)))

### JSD + PCoA ###

library(Tmisc)
library(ape)

long_jsd_pcoa <- function(data, method){
    rel_counts <- data %>% 
        filter(assembler==method) %>%
        select(species, cluster, read_count_cluster) %>% 
        group_by(cluster, species) %>% 
        summarise(count=sum(read_count_cluster)) %>% 
        ungroup() %>% 
        group_by(cluster) %>% 
        mutate(abundance=count/sum(as.numeric(count)))

    owner_info <- rel_counts %>% 
        summarize(
            cluster_owner=species[abundance==max(abundance)][1], 
            abundance=max(abundance), 
            count=max(count)) %>% 
        mutate(owner=ifelse(abundance>0.5, cluster_owner, NA)) %>%
        print()

    rel_matrix <- rel_counts %>%
        select(-count) %>%
        spread(cluster, abundance) %>%
        remove_rownames %>%
        column_to_rownames(var='species') %>%
        as.matrix()

    dist_matrix <- jsd(rel_matrix)
    res <- pcoa(dist_matrix)
    coords <- res$vectors %>% as_tibble()
    owner_info$cluster <- as.character(owner_info$cluster)

    cluster_info <- bind_cols(owner_info, coords)

    return(cluster_info)
}

long_euclidean_pcoa <- function(data, method){
    rel_counts <- data %>% 
        filter(assembler==method) %>%
        select(species, cluster, read_count_cluster) %>% 
        group_by(cluster, species) %>% 
        summarise(count=sum(read_count_cluster)) %>% 
        ungroup() %>% 
        group_by(cluster) %>% 
        mutate(abundance=count/sum(as.numeric(count)))

    owner_info <- rel_counts %>% 
        summarize(
            cluster_owner=species[abundance==max(abundance)][1], 
            abundance=max(abundance), 
            count=max(count)) %>% 
        mutate(owner=ifelse(abundance>0.5, cluster_owner, NA)) %>%
        print()

    rel_matrix <- rel_counts %>%
        select(-count) %>%
        spread(cluster, abundance) %>%
        remove_rownames %>%
        column_to_rownames(var='species') %>%
        as.matrix()

    dist_matrix <- dist(t(rel_matrix), method='euclidean')
    res <- pcoa(dist_matrix)
    coords <- res$vectors %>% as_tibble()
    owner_info$cluster <- as.character(owner_info$cluster)

    cluster_info <- bind_cols(owner_info, coords)

    return(cluster_info)
}

lsa_55_jsd <- long_jsd_pcoa(binners_count_info, 'LSA_55')
lsa_33_jsd <- long_jsd_pcoa(binners_count_info, 'LSA_33')
lsa_15_jsd <- long_jsd_pcoa(binners_count_info, 'LSA_15')

lsa_55_euclid <- long_euclidean_pcoa(binners_count_info, 'LSA_55')
lsa_33_euclid <- long_euclidean_pcoa(binners_count_info, 'LSA_33')
lsa_15_euclid <- long_euclidean_pcoa(binners_count_info, 'LSA_15')

dist_matrix_euclid <- function(data, method){
    matr <- data %>% 
        filter(assembler==method) %>%
        ungroup() %>%
        select(species, cluster, read_count_cluster) %>% 
        group_by(cluster, species) %>% 
        summarise(count=sum(read_count_cluster)) %>% 
        ungroup() %>% 
        group_by(cluster) %>% 
        mutate(abundance=count/sum(as.numeric(count))) %>%
        select(-count) %>%
        spread(cluster, abundance) %>%
        remove_rownames %>%
        column_to_rownames(var='species') %>%
        as.matrix()

    out_matrix <- dist(t(matr), method='euclidean')

    return(as.matrix(out_matrix))
}

dist_matrix_jsd <- function(data, method){
    matr <- data %>% 
        filter(assembler==method) %>%
        ungroup() %>%
        select(species, cluster, read_count_cluster) %>% 
        group_by(cluster, species) %>% 
        summarise(count=sum(read_count_cluster)) %>% 
        ungroup() %>% 
        group_by(cluster) %>% 
        mutate(abundance=count/sum(as.numeric(count))) %>%
        select(-count) %>%
        spread(cluster, abundance) %>%
        remove_rownames %>%
        column_to_rownames(var='species') %>%
        as.matrix()

    out_matrix <- jsd(matr)

    return(as.matrix(out_matrix))
}

lsa_55_euclid_matrix <- dist_matrix_euclid(binners_count_info, 'LSA_55')
lsa_33_euclid_matrix <- dist_matrix_euclid(binners_count_info, 'LSA_33')
lsa_15_euclid_matrix <- dist_matrix_euclid(binners_count_info, 'LSA_15')

lsa_55_jsd_matrix <- dist_matrix_jsd(binners_count_info, 'LSA_55')
lsa_33_jsd_matrix <- dist_matrix_jsd(binners_count_info, 'LSA_33')
lsa_15_jsd_matrix <- dist_matrix_jsd(binners_count_info, 'LSA_15')



### FIDDELING WITH STRAINS ### 
ani_tbl <- read_tsv('results/perc_ids.tab', skip=2) %>% gather(key=representative, value=ANI, -X1) %>% rename('species'=X1) %>% print() 

reps <- read_tsv('nile/data/specs/representative_genomes_specI.txt', col_names = F)

strains_tbl_low <- read_tsv('nile/data/sim/low.txt', col_names = F) %>% 
    mutate(status='representative', representative=X1, complexity='low') %>%
    rename('species'=X1)

med <- read_tsv('nile/data/sim/med.txt', col_names = F) %>% 
    mutate(status=ifelse(X1 %in% low$species, 'representative', 'strain'))

n <- 0
strains_tbl_medium <- c()
for (spec in med$X1){
  print(spec)
  if (spec %in% reps$X1){
    n <- 0
    r <- spec
    row <- c(spec, 'representative', spec)
    names(row) <- c('species', 'status', 'representative')
    strains_tbl_medium <- bind_rows(strains_tbl_medium, row)
  } else {
    row <- c(spec, paste0('strain_', n), r)
    n <- n + 1
    names(row) <- c('species', 'status', 'representative')
    strains_tbl_medium <- bind_rows(strains_tbl_medium, row)
    row <- c(r, spec)
  }
}
strains_tbl_medium <- strains_tbl_medium %>% mutate(complexity='medium')

complex <- read_tsv('nile/data/sim/complex.txt', col_names = F) %>% print()

n <- 0
strains_tbl_complex <- c()
for (spec in complex$X1){
  print(spec)
  if (spec %in% reps$X1){
    n <- 0
    r <- spec
    row <- c(spec, 'representative', spec)
    names(row) <- c('species', 'status', 'representative')
    strains_tbl_complex <- bind_rows(strains_tbl_complex, row)
  } else {
    row <- c(spec, paste0('strain_', n), r)
    n <- n + 1
    names(row) <- c('species', 'status', 'representative')
    strains_tbl_complex <- bind_rows(strains_tbl_complex, row)
    row <- c(r, spec)
  }
}
strains_tbl_complex <- strains_tbl_complex %>% mutate(complexity='complex')

strains_tbl <- bind_rows(strains_tbl_complex, strains_tbl_medium, strains_tbl_low)


taxonomy <- read_tsv('nile/data/specs/specI_taxonomy.tsv')
mapping <- read_tsv('nile/data/sim/top_70_representatives.txt', col_names=F)
mapping$X1 <- str_replace(mapping$X1, '.*\\{(.*)\\}', '\\1')

taxonomy_tbl <- inner_join(taxonomy, mapping, by=c('specI_cluster'='X1')) %>% select(X3, phylum, family)



### POOLED ASSEMBLY ANALYSIS ###

# lsa_15_500b <- read_tsv('results/new/lsa_15_500b.tsv')
# lsa_15_10kb <- read_tsv('results/new/lsa_15_10kb.tsv')
lsa_33_500b <- read_tsv('results/new/lsa_33_500b.tsv')
lsa_33_10kb <- read_tsv('results/new/lsa_33_10kb.tsv')
concoct_66_500b <- read_tsv('results/new/concoct_66_500b.tsv')
concoct_66_10kb <- read_tsv('results/new/concoct_66_10kb.tsv')

megahit_30_samples_500b <- read_tsv('results/new/megahit_30_samples_500b.tsv')
megahit_30_samples_10kb  <- read_tsv('results/new/megahit_30_samples_10kb.tsv')
megahit_10_samples_500b  <- read_tsv('results/new/megahit_10_samples_500b.tsv')
megahit_10_samples_10kb  <- read_tsv('results/new/megahit_10_samples_10kb.tsv')

lsa_33 <- bind_rows('500b'=lsa_33_500b, '10kb'=lsa_33_10kb, .id = 'cutoff')
concoct_66 <- bind_rows('500b'=concoct_66_500b, '10kb'=concoct_66_10kb, .id = 'cutoff')
# lsa_15 <- bind_rows('500b'=lsa_15_500b, '10kb'=lsa_15_10kb, .id = 'cutoff')
megahit_30 <- bind_rows('500b'=megahit_30_samples_500b, '10kb'=megahit_30_samples_10kb, .id = 'cutoff')
megahit_10 <- bind_rows('500b'=megahit_10_samples_500b, '10kb'=megahit_10_samples_10kb, .id = 'cutoff')

pooled_tbl <- bind_rows(lsa_33, concoct_66, megahit_10, megahit_30) %>% rename('method'=sample, 'total_length_10000'=total_length_1000, 'total_length_1000'=total_length_10000)

# colnames(pooled_tbl)[8:9] <- c('misassembled_contigs_length', 'misassemblies')