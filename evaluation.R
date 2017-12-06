
library(tidyverse)

setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')
source('nile/scripts/plotting/plots.R')
source('nile/scripts/plotting/data.R')


### SINGLE SAMPLE ASSEMBLY ###

comp_order <- c('low', 'medium', 'complex')
assembly_info_mc$complexity <- factor(assembly_info_mc$complexity, levels = comp_order)
assembly_info_mc <- assembly_info_mc %>% 
    arrange(desc(complexity))
assembly_info_mc <- assembly_info_mc %>% mutate(misassembled_fraction=(missassembled_contigs_length/assembly_len) * 100)

# spades_samples <- unique(assembly_info_mc$sample[assembly_info_mc$assembler=='SPADES'])
# megahit_samples <- unique(assembly_info_mc$sample[assembly_info_mc$assembler=='MEGAHIT'])
# both <- intersect(spades_samples, megahit_samples)

this <- assembly_info_mc %>% filter(cutoff == '500b', complexity == 'low')
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.5, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)
this <- assembly_info_mc %>% filter(cutoff == '500b', complexity == 'medium')
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.9, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)
this <- assembly_info_mc %>% filter(cutoff == '500b', complexity == 'complex')
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.9, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)

#####
### Different whiskers and quantiles for boxes 
#####
this <- assembly_info_mc %>% filter(cutoff == '500b')
title <- 'Genome Fraction of Assemblies for Different Complexities (500bp Cutoff)'
boxPlot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', title=title, colorscheme=material)
boxPlot(this, this$complexity, this$misassemblies, fill=this$assembler, xlab='Assembler', ylab='Misassemblies', flab='Complexity', title=title, colorscheme=material)
boxPlot(this, this$complexity, this$NGA50, fill=this$assembler, xlab='Assembler', ylab='NGA50', flab='Complexity', title=title, colorscheme=material)
boxPlot(this, this$complexity, this$LGA50, fill=this$assembler, xlab='Assembler', ylab='LGA50', flab='Complexity', title=title, colorscheme=material)


randombeePlot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', title=title, colorscheme=material)
randombeePlot(this, this$complexity, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', title=title, colorscheme=material)
randombeePlot(this, this$complexity, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', title=title, colorscheme=material)
randombeePlot(this, this$complexity, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', title=title, colorscheme=material)


this <- assembly_info_mc %>% filter(cutoff == '10kb')
title <- 'Genome Fraction of Assemblies for Different Complexities (10kp Cutoff)'
boxPlot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', title=title, colorscheme=material)
randombeePlot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', title=title, colorscheme=material)

# this <- assembly_info_mc %>% filter(cutoff == '500b')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 500bp Cutoff)'
# qqPlot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qqPlot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qqPlot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])

# this <- assembly_info_mc %>% filter(cutoff == '10kb')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 10kb Cutoff)'
# qqPlot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qqPlot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qqPlot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])


nq <- 1000
p <- (1 : nq) / nq - 0.5 / nq
this <- assembly_info_mc %>% filter(cutoff == '500b')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')
title <- 'Genome Fraction: MEGAHIT Quantiles vs metaSPAdes Quantiles (500bp Cutoff)'
qqPlotSS(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
title <- 'Misassemblies: MEGAHIT Quantiles vs metaSPAdes Quantiles (500bp Cutoff)'
qqPlotSS(x=mega, y=spades, variable='misassemblies', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
title <- 'NGA50: MEGAHIT Quantiles vs metaSPAdes Quantiles (500bp Cutoff)'
qqPlotSS(x=mega, y=spades, variable='NGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
title <- 'LGA50: MEGAHIT Quantiles vs metaSPAdes Quantiles (500bp Cutoff)'
qqPlotSS(x=mega, y=spades, variable='LGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
title <- 'Misassembled contig fraction: MEGAHIT Quantiles vs metaSPAdes Quantiles (500bp Cutoff)'
qqPlotSS(x=mega, y=spades, variable='misassembled_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)


title <- 'Genome Fraction: MEGAHIT Quantiles vs metaSPAdes Quantiles (10kb Cutoff)'
this <- assembly_info_mc %>% filter(cutoff == '10kb')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')
qqPlotSS(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', dist=p, legend_title='Complexity', legend_labels=c('all', 'low', 'medium', 'complex'), title=title, colorscheme=material)

this <- assembly_info_mc %>% filter(cutoff == '500b')
title <- 'ECDF Plot: Genome Fraction' 
ecdfPlot(this, this$genome_fraction, col=this$complexity, line=this$assembler, title=title, xlab='Cumulative Density', ylab='Genome fraction [%]', colorscheme=material)
title <- 'ECDF Plot: LGA50' 
ecdfPlot(this, this$LGA50, col=this$complexity, line=this$assembler, title=title, xlab='Cumulative Density', ylab='LGA50', colorscheme=material)
title <- 'ECDF Plot: NGA50' 
ecdfPlot(this, this$NGA50, col=this$complexity, line=this$assembler, title=title, xlab='Cumulative Density', ylab='NGA50', colorscheme=material)
title <- 'ECDF Plot: Misassemblies' 
ecdfPlot(this, this$misassemblies, col=this$complexity, line=this$assembler, title=title, xlab='Cumulative Density', ylab='Misassemblies', colorscheme=material)
title <- 'ECDF Plot: Total length (> 10kb)' 
ecdfPlot(this, this$total_length_10000bp, col=this$complexity, line=this$assembler, title=title, xlab='Cumulative Density', ylab='Total length [bp]', colorscheme=material)


### EVALIUATION OF EFFECTS OF STRAIN PRESENCE ###

as_bin_tbl <- inner_join(as_bin_tbl, strains_tbl)
as_bin_tbl$read_count[as_bin_tbl$complexity == 'low'] <- as_bin_tbl$read_count[as_bin_tbl$complexity == 'low'] * 2
as_bin_tbl$genome_coverage[as_bin_tbl$complexity == 'low'] <- as_bin_tbl$genome_coverage[as_bin_tbl$complexity == 'low'] * 2
order_vec <- c('low', 'medium', 'complex')
as_bin_tbl$complexity <- factor(as_bin_tbl$complexity, levels=order_vec)

this <- as_bin_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>0.005, cutoff=='10kb', sample=='SAMEA3136644')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', size=3)

this <- as_bin_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>10, cutoff=='10kb', sample=='SAMEA3136644', representative %in% c('357276.PRJNA232731', '1121115.PRJNA195783', '245018.PRJNA45957', '411479.PRJNA18195', '1680.PRJNA240293', '657322.PRJNA39151'))
logPlot(this, this$genome_coverage, this$genome_fraction, this$representative, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', alpha=1, size=4, slab='Complexity')


this <- as_bin_tbl %>% filter(assembler=='MEGAHIT', genome_coverage>0.005, (genome_coverage<100 | genome_fraction>50), cutoff=='500b')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=2, colorscheme=material, alpha=0.5)
this <- as_bin_tbl %>% filter(assembler=='MEGAHIT', genome_coverage>20, (genome_coverage<100 | genome_fraction>50), cutoff=='500b')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=2, colorscheme=material, alpha=0.5)

# Select intersection of samples and SS assemblers
as_bin_tbl <- as_bin_tbl %>% filter(assembler %in% c('MEGAHIT', 'SPADES'), cutoff=='500b')
complex_sample <- unique(as_bin_tbl$sample[as_bin_tbl$complexity=='complex'])
medium_sample <- unique(as_bin_tbl$sample[as_bin_tbl$complexity=='medium'])
low_sample <- unique(as_bin_tbl$sample[as_bin_tbl$complexity=='low'])
inter_sample <- intersect(low_sample,medium_sample)
inter_tbl <- as_bin_tbl %>% filter(sample %in% inter_sample)

# Select species in upper range but below 85% genome fraction as outliers
outlier_species <- unique(inter_tbl %>% filter(genome_fraction<85, genome_coverage>20, complexity=='low') %>% select(representative))
robust_species <- unique(inter_tbl %>% filter(complexity=='low', !(representative %in% outlier_species$representative)) %>% select(representative)) 
robust_tbl <- inter_tbl %>% filter(genome_coverage>20, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', (representative %in% robust_species$representative))
this <- robust_tbl %>% filter(complexity!='complex')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=2, colorscheme=material, alpha=0.5)


# Calculate loss of coverage


robust_tbl <- mutate(LOGF=)



### BINNERS AND SINGE SAMPLE COMPARISON ###

comp_order <- c('low', 'medium', 'complex')
as_bin_tbl$complexity <- factor(as_bin_tbl$complexity, levels = comp_order)
as_bin_tbl <- as_bin_tbl %>% arrange(desc(complexity))

as_bin_tbl <- as_bin_tbl %>% mutate(abundance_bins=cut(rel_abundance, breaks=c(-Inf, 0.00001, 0.0001, 0.001, 0.01, 0.1, +Inf), labels=c('< 0.00001', '0.00001 - 0.0001', '0.0001 - 0.001', '0.001 - 0.01', '0.01 - 0.1', '> 0.1'))) 

this <- as_bin_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '500b')
title <- 'Comparison of Different Assembly Strategies (500bp Cutoff)'
boxPlot(this, this$assembler, this$genome_fraction, fill=this$abundance_bins, xlab='Assembler', ylab='Genome fraction [%]', flab='Relative Abundance', title=title, colorscheme=colorful)
boxPlot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Relative Abundance', ylab='Genome fraction [%]', flab='Assembler', title=title, colorscheme=colorful)
boxPlot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', title=title, colorscheme=colorful)
barPlotDodge(this, this$assembler, this$NGA50, fill=this$assembler, xlab='Assembler', ylab='NGA50', flab='Complexity', title=title, colorscheme=colorful)
barPlotDodge(this, this$assembler, this$LGA50, fill=this$assembler, xlab='Assembler', ylab='LGA50', flab='Complexity', title=title, colorscheme=colorful)

this <- as_bin_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '10kb')
title <- 'Comparison of Different Assembly Strategies (10kb Cutoff)'
boxPlot(this, this$assembler, this$genome_fraction, fill=this$abundance_bins, xlab='Assembler', ylab='Genome fraction [%]', flab='Relative Abundance', title=title, colorscheme=colorful)
boxPlot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Relative Abundance', ylab='Genome fraction [%]', flab='Assembler', title=title, colorscheme=colorful)
boxPlot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', title=title, colorscheme=colorful)
barPlotDodge(this, this$assembler, this$NGA50, fill=this$assembler, xlab='Assembler', ylab='NGA50', flab='Complexity', title=title, colorscheme=colorful)
barPlotDodge(this, this$assembler, this$LGA50, fill=this$assembler, xlab='Assembler', ylab='LGA50', flab='Complexity', title=title, colorscheme=colorful)



title <- 'ECDF Plot: Genome Fraction' 
this <- as_bin_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '500b')
ecdfPlot(this, this$genome_fraction, col=this$assembler, title=title, colorscheme=colorful)

### EVALUATION OF BINNERS/CLUSTERING ###

title <- 'CONCOCT: Number of reads per Cluster for each Species'
this <- concoct_counts_medium %>% filter(sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
barPlotCount(this, this$cluster, this$read_count, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many)
barPlotCount(this, this$cluster, this$read_abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many)

title <- 'LSA: Number of reads per Cluster for each Species'
this <- lsa_counts_medium %>% filter(sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
barPlotCount(this, this$cluster, this$read_count, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many)
barPlotCount(this, this$cluster, this$read_abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many)



