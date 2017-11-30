
library(tidyverse)

setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')
source('nile/scripts/plotting/plots.R')
source('nile/scripts/plotting/data.R')


comp_order <- c('low', 'medium', 'complex')
assembly_info_mc$complexity <- factor(assembly_info_mc$complexity, levels = comp_order)
assembly_info_mc <- assembly_info_mc %>% 
    arrange(desc(complexity))


this <- assembly_info_mc %>% filter(cutoff == '500b')
boxPlot2(this, this$assembler, this$genome_fraction, fill=this$complexity, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity')
boxPlot2(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity')

this <- assembly_info_mc %>% filter(cutoff == '10kb')
boxPlot2(this, this$assembler, this$genome_fraction, fill=this$complexity, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity')
boxPlot2(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity')

this <- assembly_info_mc %>% filter(cutoff == '500b')
qqPlot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity')
qqPlot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity')
qqPlot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity')

nq <- 1000
p <- (1 : nq) / nq - 0.5 / nq
this <- assembly_info_mc %>% filter(cutoff == '500b')
that <- this
that$complexity <- 'all'
this <- this %>% bind_rows(that)
this$complexity <- factor(this$complexity, levels=c('all', 'low', 'medium', 'complex'))
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')
qqPlotSS(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('all', 'low', 'medium', 'complex'))

ecdfPlot(this, this$genome_fraction, col=this$assembler)
ecdfPlot(this, this$genome_fraction, col=this$complexity, line=this$assembler)
ecdfPlot(this, log10(this$LGA50), col=this$complexity, line=this$assembler)
ecdfPlot(this, log10(this$NGA50), col=this$complexity, line=this$assembler)
ecdfPlot(this, log10(this$misassemblies), col=this$complexity, line=this$assembler)
ecdfPlot(this, this$total_length_10000bp, col=this$complexity, line=this$assembler)


### BINNERS AND SINGE SAMPLE COMPARISON ###

comp_order <- c('low', 'medium', 'complex')
as_bin_tbl$complexity <- factor(as_bin_tbl$complexity, levels = comp_order)
as_bin_tbl <- as_bin_tbl %>% arrange(desc(complexity))

this <- as_bin_tbl %>% mutate(abundance_bins=cut(rel_abundance, breaks=c(-Inf, 0.00001, 0.0001, 0.001, 0.01, 0.1, +Inf), labels=c('< 0.00001', '0.00001 - 0.0001', '0.0001 - 0.001', '0.001 - 0.01', '0.01 - 0.1', '> 0.1'))) 
boxPlot2(this, this$assembler, this$genome_fraction, fill=this$abundance_bins, xlab='Assembler', ylab='Genome fraction [%]', flab='Relative Abundance')
boxPlot2(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Relative Abundance', ylab='Genome fraction [%]', flab='Assembler')

this <- as_bin_tbl %>% filter(complexity == 'medium')
boxPlot2(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity')

ecdfPlot(this, this$genome_fraction, col=this$assembler)
ecdfPlot(this, log10(this$LGA50), col=this$assembler)
ecdfPlot(this, log10(this$NGA50), col=this$assembler)
ecdfPlot(this, this$total_length_10000bp, col=this$assembler)


### EVALUATION OF BINNERS/CLUSTERING ###
this <- concoct_counts_medium %>% filter(sample %in% c('SAMEA3136631'))
barPlotCount(this, this$cluster, this$read_count, fill=this$species, xlab='Cluster', ylab='Read Abundance')
barPlotCount(this, this$cluster, this$read_abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance')

this <- lsa_counts_medium %>% filter(sample %in% c('SAMEA3136631'))
barPlotCount(this, this$cluster, this$read_count, fill=this$species, xlab='Cluster', ylab='Read Abundance')
barPlotCount(this, this$cluster, this$read_abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance')

this <- binners_count_info %>% filter(assembler=='CONCOCT', rel_abundance>0, max_rel_abundance>0)
loglogPlot(this, this$rel_abundance, this$max_rel_abundance, this$sample, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Sample')
this <- binners_count_info %>% filter(assembler=='LSA', rel_abundance>0, max_rel_abundance>0)
loglogPlot(this, this$rel_abundance, this$max_rel_abundance, this$sample, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Sample')


### EVALIUATION OF EFFECTS OF STRAIN PRESENCE ###
as_bin_tbl <- inner_join(as_bin_tbl, strains_tbl)
as_bin_tbl$read_count[as_bin_tbl$complexity == 'low'] <- as_bin_tbl$read_count[as_bin_tbl$complexity == 'low'] * 2
as_bin_tbl$genome_coverage[as_bin_tbl$complexity == 'low'] <- as_bin_tbl$genome_coverage[as_bin_tbl$complexity == 'low'] * 2

this <- as_bin_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>0.005, cutoff=='10kb', sample=='SAMEA3136644')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction', clab='Species', size=3)

this <- as_bin_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>10, cutoff=='10kb', sample=='SAMEA3136644', representative %in% c('357276.PRJNA232731', '1121115.PRJNA195783', '245018.PRJNA45957', '411479.PRJNA18195', '1680.PRJNA240293', '657322.PRJNA39151'))
logPlot(this, this$genome_coverage, this$genome_fraction, this$representative, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction', clab='Species', alpha=1, size=4, slab='Complexity')
logPlot(this, this$rel_abundance, this$genome_fraction, this$representative, this$complexity, xlab='Relative Abundance', ylab='Genome Fraction', clab='Species', alpha=1, size=4, slab='Complexity')
logPlot(this, this$read_count, this$genome_fraction, this$representative, this$complexity, xlab='Read Count', ylab='Genome Fraction', clab='Species', alpha=1, size=4, slab='Complexity')
