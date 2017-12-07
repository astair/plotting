
library(tidyverse)


setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')
source('nile/scripts/plotting/plots.R')
source('nile/scripts/plotting/data.R')

print(the_total_tbl)

unique(the_total_tbl$assembler)

### SINGLE SAMPLE ASSEMBLY ###

ass_order <- c('MEGAHIT', 'SPADES', 'CONCOCT', 'LSA_33', 'LSA_55')
the_total_tbl$assembler <- factor(the_total_tbl$assembler, levels = ass_order)
comp_order <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels = comp_order)
the_total_tbl <- the_total_tbl %>% 
    arrange(desc(complexity))
the_total_tbl <- the_total_tbl %>% mutate(misassembled_fraction=(missassembled_contigs_length/assembly_len))

# spades_samples <- unique(the_total_tbl$sample[the_total_tbl$assembler=='SPADES'])
# megahit_samples <- unique(the_total_tbl$sample[the_total_tbl$assembler=='MEGAHIT'])
# both <- intersect(spades_samples, megahit_samples)

this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'low', assembler %in% c('MEGAHIT', 'SPADES'))

# pdf("GFvsMISS.pdf", width=6, height=4)
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.5, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)
this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'medium', assembler %in% c('MEGAHIT', 'SPADES'))
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.9, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)
this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'complex', assembler %in% c('MEGAHIT', 'SPADES'))
dotPlotSmooth(this, this$misassembled_fraction, this$genome_fraction, this$assembler, span=0.9, alpha=0.2, xlab='Fraction misassembled contig length [%]', ylab='Genome fraction [%]', size=1, colorscheme=material)
# dev.off()

#####
### Different whiskers and quantiles for boxes 
#####

this <- the_total_tbl %>% filter(cutoff == '500b', assembler %in% c('MEGAHIT', 'SPADES'))

# pdf("ss_assembly_box_500bp.pdf", width=6, height=4)
boxPlot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
boxPlot(this, this$complexity, this$misassemblies, fill=this$assembler, xlab='Assembler', ylab='Misassemblies', flab='Complexity', colorscheme=material)
boxPlot(this, this$complexity, this$NGA50, fill=this$assembler, xlab='Assembler', ylab='NGA50', flab='Complexity', colorscheme=material)
boxPlot(this, this$complexity, this$LGA50, fill=this$assembler, xlab='Assembler', ylab='LGA50', flab='Complexity', colorscheme=material)
# dev.off()

# pdf("ss_assembly_bee_500bp.pdf", width=6, height=4)
randombeePlot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material)
randombeePlot(this, this$complexity, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', colorscheme=material)
randombeePlot(this, this$complexity, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', colorscheme=material)
randombeePlot(this, this$complexity, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', colorscheme=material)
# dev.off()

# pdf("ss_assembly_10kb.pdf", width=6, height=4)
this <- the_total_tbl %>% filter(cutoff == '10kb', assembler %in% c('MEGAHIT', 'SPADES'))
boxPlot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombeePlot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material)
# dev.off()

# this <- the_total_tbl %>% filter(cutoff == '500b')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 500bp Cutoff)'
# qqPlot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qqPlot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qqPlot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])

# this <- the_total_tbl %>% filter(cutoff == '10kb')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 10kb Cutoff)'
# qqPlot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qqPlot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qqPlot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])


nq <- 1000
p <- (1 : nq) / nq - 0.5 / nq
this <- the_total_tbl %>% filter(cutoff == '500b')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')

# pdf("ss_assembly_QQ_500bp.pdf", width=6, height=4)

title <- 'Genome Fraction'
qqPlotSS(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'Misassemblies'
qqPlotSS(x=mega, y=spades, variable='misassemblies', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'NGA50'
qqPlotSS(x=mega, y=spades, variable='NGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'LGA50'
qqPlotSS(x=mega, y=spades, variable='LGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'Misassembled contig fraction'
qqPlotSS(x=mega, y=spades, variable='misassembled_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
# dev.off()

title <- 'Genome Fraction'
this <- the_total_tbl %>% filter(cutoff == '10kb')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')

# pdf("ss_assembly_QQ_10kb.pdf", width=6, height=4)
qqPlotSS(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('all', 'low', 'medium', 'complex'), title=title, colorscheme=material)
# dev.off()


this <- the_total_tbl %>% filter(cutoff == '500b', assembler %in% c('MEGAHIT', 'SPADES'))

# pdf("ss_assembly_ECDF_500bp.pdf", width=6, height=4)

ecdfPlot(this, this$genome_fraction, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Genome fraction [%]', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: LGA50' 
ecdfPlot(this, this$LGA50, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='LGA50', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: NGA50' 
ecdfPlot(this, this$NGA50, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='NGA50', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: Misassemblies' 
ecdfPlot(this, this$misassemblies, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Misassemblies', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: Total length (> 10kb)' 
ecdfPlot(this, this$total_length_10000bp, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Total length [bp]', llab='Assembler', clab='Complexity', colorscheme=material)
# dev.off()


### EVALIUATION OF EFFECTS OF STRAIN PRESENCE ###

the_total_tbl <- inner_join(the_total_tbl, strains_tbl)
the_total_tbl$read_count[the_total_tbl$complexity == 'low'] <- the_total_tbl$read_count[the_total_tbl$complexity == 'low'] * 2
the_total_tbl$genome_coverage[the_total_tbl$complexity == 'low'] <- the_total_tbl$genome_coverage[the_total_tbl$complexity == 'low'] * 2
order_vec <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels=order_vec)

# this <- the_total_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>0.005, cutoff=='10kb', sample=='SAMEA3136644')
# logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', size=3)

# this <- the_total_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>10, cutoff=='10kb', sample=='SAMEA3136644', representative %in% c('357276.PRJNA232731', '1121115.PRJNA195783', '245018.PRJNA45957', '411479.PRJNA18195', '1680.PRJNA240293', '657322.PRJNA39151'))
# logPlot(this, this$genome_coverage, this$genome_fraction, this$representative, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', alpha=1, size=4, slab='Complexity')

# pdf("ss_assembly_strains.pdf", width=6, height=4)

this <- the_total_tbl %>% filter(assembler=='MEGAHIT', genome_coverage>0.005, (genome_coverage<100 | genome_fraction>50), cutoff=='500b')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=material, alpha=0.5)
this <- the_total_tbl %>% filter(assembler=='MEGAHIT', genome_coverage>20, (genome_coverage<100 | genome_fraction>50), cutoff=='500b')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=material, alpha=0.5)

# Select intersection of samples and SS assemblers
the_total_tbl <- the_total_tbl %>% filter(assembler %in% c('MEGAHIT', 'SPADES'), cutoff=='500b')
complex_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='complex'])
medium_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='medium'])
low_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='low'])
inter_sample <- intersect(low_sample,medium_sample)
inter_tbl <- the_total_tbl %>% filter(sample %in% inter_sample)

# Select species in upper range but below 85% genome fraction as outliers
outlier_species <- unique(inter_tbl %>% filter(genome_fraction<85, genome_coverage>20, complexity=='low') %>% select(representative))
robust_species <- unique(inter_tbl %>% filter(complexity=='low', !(representative %in% outlier_species$representative)) %>% select(representative)) 
robust_tbl <- inter_tbl %>% filter(genome_coverage>20, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', (representative %in% robust_species$representative))
this <- robust_tbl %>% filter(complexity!='complex')
logPlot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=material, alpha=0.5)


# Calculate loss of coverage
rep_coverage <- robust_tbl %>% 
    filter(complexity=='low') %>% 
    select(representative, sample, genome_fraction) %>% 
    rename('reference_coverage'=genome_fraction)
loss_tbl <- inner_join(rep_coverage, robust_tbl) %>% 
    mutate(coverage_loss=reference_coverage-genome_fraction) %>%
    inner_join(ani_tbl)

this <- loss_tbl %>% filter(ANI > 0.9)
dotPlot(this, this$ANI, this$coverage_loss, this$assembler, xlab='ANI [%]', ylab='Loss of Genome Fraction [%]', clab='Assembler', size=1, colorscheme=material, alpha=0.5)
dotPlotSmooth(this, this$ANI, this$coverage_loss, this$assembler, xlab='ANI [%]', ylab='Loss of Genome Fraction [%]', clab='Assembler', size=1, colorscheme=material, alpha=0.5)

# dev.off()


### BINNERS AND SINGE SAMPLE COMPARISON ###

comp_order <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels = comp_order)
the_total_tbl <- the_total_tbl %>% arrange(desc(complexity))

the_total_tbl <- the_total_tbl %>% mutate(abundance_bins=cut(genome_coverage, breaks=c(-Inf, 1, 2.5, 5, 7.5, 10, +Inf), labels=c('< 1x', '1x - 2.5x', '2.5x - 5x', '5x - 7.5x', '7.5x - 10x', '> 10x'))) 

# pdf("binnes_vs_ss_assembly_500bp.pdf", width=6, height=4)

this <- the_total_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '500b')

boxPlot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Fold Genome Coverage', ylab='Genome fraction [%]', flab='Assembler', colorscheme=material)
boxPlot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombeePlot(this, this$assembler, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material)
randombeePlot(this, this$assembler, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', colorscheme=material)
randombeePlot(this, this$assembler, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', colorscheme=material)
randombeePlot(this, this$assembler, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', colorscheme=material)

# dev.off()

# pdf("binnes_vs_ss_assembly_10kb.pdf", width=6, height=4)

this <- the_total_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '10kb')

boxPlot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Fold Genome Coverage', ylab='Genome fraction [%]', flab='Assembler', colorscheme=material)
boxPlot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombeePlot(this, this$assembler, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material)
randombeePlot(this, this$assembler, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', colorscheme=material)
randombeePlot(this, this$assembler, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', colorscheme=material)
randombeePlot(this, this$complexity, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', colorscheme=material)

# dev.off()

title <- 'ECDF Plot: Genome Fraction' 
this <- the_total_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '500b')
ecdfPlot(this, this$genome_fraction, col=this$assembler, title=title, colorscheme=material)

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



