
library(tidyverse)


setwd('~/Documents/Msc-Biotechnologie/masterarbeit-zeller/')
source('nile/scripts/plotting/plots.R')
source('nile/scripts/plotting/data_new.R')

# the_total_tbl[,1:15]

# unique(the_total_tbl$assembler)

### SAMPLES AND STRAINS ###

#     group_by(cluster) %>% 
#     mutate(cluster_sum = sum(read_count)) %>%
#     arrange(desc(cluster_sum))
# this$cluster <- factor(this$cluster, levels = unique(this$cluster))

used_samples <- unique(the_total_tbl$sample)

info <- inner_join(strains_tbl, info_all) %>% inner_join(taxonomy_tbl, by=c('representative'='X3')) %>% filter(sample %in% used_samples)
info$phylum[info$phylum=='NA Bacteroidales fam. [C Bacteroidaceae/Porphyromonadaceae]'] <- 'Bacteroidales'
info$phylum[info$phylum=='NA Bacteria fam. incertae sedis'] <- 'Incertae sedis'
info$phylum[info$phylum=='NA Bacteria phylum incertae sedis'] <- 'Incertae sedis'
info$phylum <- str_replace(info$phylum, '\\d+ (\\w+)', '\\1')
status_order <- c('representative', 'strain_0', 'strain_1', 'strain_2', 'strain_3', 'strain_4', 'strain_5', 'strain_6')
info$status <- factor(info$status, levels=status_order)

# pdf("datasets_overview.pdf", width=12, height=3)

this <- info %>% 
    filter(complexity=='low') 
spec_order <- this %>% group_by(representative) %>% summarise(median_abundance=median(rel_abundance)) %>% arrange(desc(median_abundance))

colors <- c('#Ae0031', '#2d62a3', '#ffbb00', '#198749', '#6c3483', '#Ff7708', '#15889C')

this$representative <- factor(this$representative, levels=spec_order$representative)
this$species <- factor(this$species, levels=unique(this$species))
this$phylum <- factor(this$phylum, levels=unique(this$phylum))

this <- this %>% filter(representative %in% c('420247.PRJNA18653', '511680.PRJNA28999', '592028.PRJNA33143', '702459.PRJNA42863'))

logbox_plot(this, this$representative, this$rel_abundance, this$phylum, this$status, out_size=0.5, ylim=c(NA, 1), colorscheme=colors, xlab='Species', ylab='Relative Abundance', flab='Phylum', blank_x=TRUE)



colors <- c('#Ae0031', '#2d62a3', '#198749', '#ffbb00', '#Ff7708', '#15889C')

this$representative <- factor(this$representative, levels=unique(this$representative))
this$species <- factor(this$species, levels=unique(this$species))
this$phylum <- factor(this$phylum, levels=unique(this$phylum))
logbox_plot(this, this$representative, this$rel_abundance, this$phylum, this$status, out_size=0.5, ylim=c(NA, 1), colorscheme=colors, xlab='Species', ylab='Relative Abundance', flab='Phylum', blank_x=TRUE)

this <- info %>% 
    filter(complexity=='complex') %>% 
    group_by(representative) %>%
    mutate(mean_abundance=median(rel_abundance)) %>%
    arrange(desc(mean_abundance))

unique(this$species)

colors <- c('#2d62a3', '#198749', '#Ae0031', '#ffbb00', '#Ff7708', '#6c3483', '#15889C')

this$representative <- factor(this$representative, levels=unique(this$representative))
this$species <- factor(this$species, levels=unique(this$species))
this$phylum <- factor(this$phylum, levels=unique(this$phylum))
logbox_plot(this, this$representative, this$rel_abundance, this$phylum, this$status, out_size=0.5, ylim=c(NA, 1), colorscheme=colors, xlab='Species', ylab='Relative Abundance', flab='Phylum', blank_x=TRUE)

# dev.off()
min(info$rel_abundance[info$rel_abundance!=0])

length(unique(info$family))



### SINGLE SAMPLE ASSEMBLY ###

ass_order <- c('MEGAHIT', 'SPADES', 'CONCOCT', 'LSA_15', 'LSA_33', 'LSA_55')
the_total_tbl$assembler <- factor(the_total_tbl$assembler, levels = ass_order)
comp_order <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels = comp_order)
the_total_tbl <- the_total_tbl %>% 
    arrange(desc(complexity))
the_total_tbl <- the_total_tbl %>% 
    mutate(misassembled_genome_fraction=(misassembled_contig_length/assembly_len))
the_total_tbl <- the_total_tbl %>% 
    mutate(misassembled_fraction=(misassembled_contig_length/total_length))
the_total_tbl <- the_total_tbl %>% 
    mutate(aligned_fraction=(total_aligned_length/total_length))
the_total_tbl <- the_total_tbl %>% 
    mutate(mis_fraction=(misassemblies/total_length))


# spades_samples <- unique(the_total_tbl$sample[the_total_tbl$assembler=='SPADES'])
# megahit_samples <- unique(the_total_tbl$sample[the_total_tbl$assembler=='MEGAHIT'])
# both <- intersect(spades_samples, megahit_samples)



the_total_tbl <- the_total_tbl %>% mutate(miss_gen_bins=cut(misassembled_genome_fraction, breaks=c(0, 0.1, 0.2, 0.3, 0.4), labels=c('0-10', '10-20', '20-30', '30-40'))) 
the_total_tbl <- the_total_tbl %>% mutate(cov_bins=cut(genome_coverage, breaks=c(-Inf, 10, 50, 200, +Inf), labels=c('0-10', '10-50', '50-200', '>200'))) 

# pdf("ss_misassemblies_box_500bp_new.pdf", width=3.5, height=3)

this <- the_total_tbl %>% filter(cov_bins != '0-10', assembler %in% c('MEGAHIT', 'SPADES'))
max(this$misassembled_fraction[!this$misassembled_fraction%in%c(NaN, NA)]) 

this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'low', assembler %in% c('MEGAHIT', 'SPADES'))

box_plot(this, this$cov_bins, this$misassembled_fraction, this$assembler, xlab='Genome fraction [%]', ylab='Misassembled genome fraction [%]', flab='Assembler', colorscheme=material, title='Low')
box_plot(this, this$cov_bins, this$mis_fraction, this$assembler, xlab='Genome fraction [%]', ylab='Misassembly frequency', flab='Assembler', colorscheme=material, title='Low', ylim=c(0,0.000025))

# this <- this %>% filter(genome_coverage>1)
# log_plot(this, this$genome_coverage, this$misassembled_fraction, this$assembler)
# log_plot(this, this$genome_coverage, this$misassembled_genome_fraction, this$assembler)
# log_plot(this, this$genome_fraction, this$misassembled_genome_fraction, this$assembler)
# randombee_plot(this, this$miss_gen_bins, this$genome_fraction, this$assembler, box=TRUE, xlab='Misassembled genome fraction [%]', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material, title='Low', alpha=0.4)


this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'medium', assembler %in% c('MEGAHIT', 'SPADES'))

box_plot(this, this$cov_bins, this$misassembled_fraction, this$assembler, xlab='Genome fraction [%]', ylab='Misassembled genome fraction [%]', flab='Assembler', colorscheme=material, title='Medium')
# randombee_plot(this, this$miss_gen_bins, this$genome_fraction, this$assembler, box=TRUE, xlab='Misassembled genome fraction [%]', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material, title='Medium', alpha=0.4)


this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'complex', assembler %in% c('MEGAHIT', 'SPADES'))

box_plot(this, this$cov_bins, this$misassembled_fraction, this$assembler, xlab='Genome fraction [%]', ylab='Misassembled genome fraction [%]', flab='Assembler', colorscheme=material, title='Complex')
# randombee_plot(this, this$miss_gen_bins, this$genome_fraction, this$assembler, box=TRUE, xlab='Misassembled genome fraction [%]', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material, title='Complex', alpha=0.4)
        
# dev.off()


this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'medium', assembler %in% c('MEGAHIT', 'SPADES'))
dot_plot_smooth(this, this$misassembled_genome_fraction, this$genome_fraction, this$assembler, span=0.6, alpha=0.2, xlab='Misassembled genome fraction [%]', ylab='Genome fraction [%]', clab='Assembler', size=1, colorscheme=material, title='Medium')

this <- the_total_tbl %>% filter(cutoff == '500b', complexity == 'complex', assembler %in% c('MEGAHIT', 'SPADES'))
dot_plot_smooth(this, this$misassembled_genome_fraction, this$genome_fraction, this$assembler, span=0.6, alpha=0.2, xlab='Misassembled genome fraction [%]', ylab='Genome fraction [%]', clab='Assembler', size=1, colorscheme=material, title='Complex')

# dev.off()


# pdf("ss_box_genomefrac_new.pdf", width=4, height=3)

this <- the_total_tbl %>% filter(cutoff == '500b', assembler %in% c('MEGAHIT', 'SPADES'))
# box_plot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
# box_plot(this, this$complexity, this$misassemblies, fill=this$assembler, xlab='Assembler', ylab='Misassemblies', flab='Complexity', colorscheme=material)
# box_plot(this, this$complexity, this$NGA50, fill=this$assembler, xlab='Assembler', ylab='NGA50', flab='Complexity', colorscheme=material)
# box_plot(this, this$complexity, this$LGA50, fill=this$assembler, xlab='Assembler', ylab='LGA50', flab='Complexity', colorscheme=material)

# violin_plot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
# violin_plot(this, this$complexity, this$misassemblies, fill=this$assembler, xlab='Assembler', ylab='Misassemblies', flab='Complexity', colorscheme=material)

randombee_plot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.03, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material, median=T) 
# randombee_plot(this, this$complexity, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', colorscheme=material)
# randombee_plot(this, this$complexity, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', colorscheme=material)
# randombee_plot(this, this$complexity, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', colorscheme=material)
# dev.off()

# pdf("ss_assembly_various_10kb_new.pdf", width=8, height=4)
this <- the_total_tbl %>% filter(cutoff == '10kb', assembler %in% c('MEGAHIT', 'SPADES'))
# box_plot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
# violin_plot(this, this$complexity, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombee_plot(this, this$complexity, this$genome_fraction, color=this$assembler, alpha=0.03, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material, median=T)
# dev.off()

# this <- the_total_tbl %>% filter(cutoff == '500b')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 500bp Cutoff)'
# qq_plot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qq_plot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qq_plot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])

# this <- the_total_tbl %>% filter(cutoff == '10kb')
# title <- 'Genome Fraction Quantiles vs theoretical Quantiles (Normal, 10kb Cutoff)'
# qq_plot(this, this$genome_fraction, col=this$assembler, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=vibrant)
# qq_plot(this, this$genome_fraction, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])
# qq_plot(this, this$genome_fraction, shape=this$assembler, col=this$complexity, xlab='Theoretical', ylab='Genome fraction [%]', clab='Complexity', title=title, colorscheme=greentored[c(1,2,4)])

# pdf("ss_genomefrac_QQ.pdf", width=4.5, height=3)

nq <- 1000
p <- (1 : nq) / nq - 0.5 / nq
this <- the_total_tbl %>% filter(cutoff == '500b')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')

title <- 'Genome Fraction'
qq_plot_ss(x=mega, y=spades, variable='genome_fraction', split='complexity', xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=three)


nq <- 1000
p <- (1 : nq) / nq - 0.5 / nq
this <- the_total_tbl %>% filter(cutoff == '10kb')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')

title <- 'Genome Fraction'
qq_plot_ss(x=mega, y=spades, variable='genome_fraction', split='complexity', xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=three)


# dev.off()





title <- 'Misassemblies'
qq_plot_ss(x=mega, y=spades, variable='misassemblies', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'NGA50'
qq_plot_ss(x=mega, y=spades, variable='NGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

title <- 'LGA50'
qq_plot_ss(x=mega, y=spades, variable='LGA50', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)

# title <- 'Misassembled contig fraction'
# qq_plot_ss(x=mega, y=spades, variable='misassembled_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('low', 'medium', 'complex'), title=title, colorscheme=material)
# dev.off()

title <- 'Genome Fraction'
this <- the_total_tbl %>% filter(cutoff == '10kb')
order_vec <- c('low', 'medium', 'complex')
this$complexity <- factor(this$complexity, levels=order_vec)
mega <- this %>% filter(assembler == 'MEGAHIT')
spades <- this %>% filter(assembler == 'SPADES')

# pdf("ss_assembly_QQ_10kb.pdf", width=6, height=4)
qq_plot_ss(x=mega, y=spades, variable='genome_fraction', split='complexity', colors=magma, xlab='Megahit quantiles', ylab='Spades quantiles', clab='Complexity', dist=p, legend_title='Complexity', legend_labels=c('all', 'low', 'medium', 'complex'), title=title, colorscheme=material)
# dev.off()


this <- the_total_tbl %>% filter(cutoff == '500b', assembler %in% c('MEGAHIT', 'SPADES'))

# pdf("ss_assembly_ECDF_500bp.pdf", width=6, height=4)

ecdf_plot(this, this$genome_fraction, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Genome fraction [%]', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: LGA50' 
ecdf_plot(this, this$LGA50, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='LGA50', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: NGA50' 
ecdf_plot(this, this$NGA50, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='NGA50', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: Misassemblies' 
ecdf_plot(this, this$misassemblies, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Misassemblies', llab='Assembler', clab='Complexity', colorscheme=material)
title <- 'ECDF Plot: Total length (> 10kb)' 
ecdf_plot(this, this$total_length_10000bp, col=this$complexity, line=this$assembler, xlab='Cumulative Density', ylab='Total length [bp]', llab='Assembler', clab='Complexity', colorscheme=material)
# dev.off()


### EVALIUATION OF EFFECTS OF STRAIN PRESENCE ###

the_total_tbl <- inner_join(the_total_tbl, strains_tbl)
the_total_tbl$read_count[the_total_tbl$complexity == 'low'] <- the_total_tbl$read_count[the_total_tbl$complexity == 'low'] * 2
the_total_tbl$genome_coverage[the_total_tbl$complexity == 'low'] <- the_total_tbl$genome_coverage[the_total_tbl$complexity == 'low'] * 2
order_vec <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels=order_vec)

# this <- the_total_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>0.005, cutoff=='10kb', sample=='SAMEA3136644')
# log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', size=3)

# this <- the_total_tbl %>% filter(complexity %in% c('medium', 'low'), assembler=='MEGAHIT', genome_coverage>10, cutoff=='10kb', sample=='SAMEA3136644', representative %in% c('357276.PRJNA232731', '1121115.PRJNA195783', '245018.PRJNA45957', '411479.PRJNA18195', '1680.PRJNA240293', '657322.PRJNA39151'))
# log_plot(this, this$genome_coverage, this$genome_fraction, this$representative, this$complexity, xlab='Genome Coverage', ylab='Genome Fraction [%]', clab='Species', alpha=1, size=4, slab='Complexity')

# pdf("ss_log_new.pdf", width=4.5, height=3)

this <- the_total_tbl %>% filter(genome_coverage>0.05, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', complexity=='low')
log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=three[1], alpha=0.3)
this <- the_total_tbl %>% filter(genome_coverage>0.05, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', complexity=='medium')
log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=three[2], alpha=0.3)
this <- the_total_tbl %>% filter(genome_coverage>0.05, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', complexity=='complex')
log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=three[3], alpha=0.3)

# dev.off()

# this <- the_total_tbl %>% filter(genome_coverage>0.005, (genome_coverage<100 | genome_fraction>50), cutoff=='500b', complexity=='medium')
# loghex_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity')

# this <- the_total_tbl %>% filter(assembler=='MEGAHIT', genome_coverage>20, (genome_coverage<100 | genome_fraction>50), cutoff=='500b')
# log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=material, alpha=0.3)

# Select intersection of samples and SS assemblers
the_total_tbl <- the_total_tbl %>% filter(assembler %in% c('MEGAHIT', 'SPADES'), cutoff=='500b')
complex_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='complex'])
medium_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='medium'])
low_sample <- unique(the_total_tbl$sample[the_total_tbl$complexity=='low'])

spades_sample <- unique(the_total_tbl$sample[the_total_tbl$assembler=='SPADES'])
mega_sample <- unique(the_total_tbl$sample[the_total_tbl$assembler=='MEGAHIT'])

inter_sample <- intersect(low_sample,medium_sample) %>% intersect(spades_sample) %>% intersect(mega_sample)
inter_tbl <- the_total_tbl %>% filter(sample %in% inter_sample)

# Select species in upper range but below 85% genome fraction as outliers
outlier_species <- unique(inter_tbl %>% filter(genome_fraction<85, genome_coverage>10, complexity=='low') %>% select(representative))
robust_species <- unique(inter_tbl %>% filter(complexity=='low', !(representative %in% outlier_species$representative)) %>% select(representative)) 
robust_tbl <- inter_tbl %>% filter(genome_coverage>10, cutoff=='500b', (representative %in% robust_species$representative))

pdf("strain_loss_robust_new.pdf", width=3, height=3)
this <- robust_tbl 
log_plot(this, this$genome_coverage, this$genome_fraction, this$complexity, xlab='Fold Genome Coverage', ylab='Genome Fraction [%]', clab='Complexity', size=1, colorscheme=three, alpha=0.5)
dev.off()

# Calculate loss of coverage
rep_coverage <- robust_tbl %>% 
    filter(complexity=='low') %>% 
    select(representative, sample, genome_fraction) %>% 
    rename('reference_coverage'=genome_fraction)
loss_tbl <- inner_join(rep_coverage, robust_tbl) %>% 
    mutate(coverage_loss=reference_coverage-genome_fraction) %>%
    inner_join(ani_tbl)
loss_tbl$ANI[is.na(loss_tbl$ANI)] <- 1
loss_tbl <- loss_tbl %>% mutate(abundance_bins=cut(genome_coverage, breaks=c(10, 100, +Inf), labels=c('10x - 100x', '> 100x'))) 


# box_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, xlab='Fold Genome Coverage [%]', ylab='Loss of Genome Fraction [%]', flab='Assembler', colorscheme=material)

# Single Reps
single <- c()
for (rep in unique(loss_tbl$representative)){
    if ( length( unique( (loss_tbl%>%filter(representative==rep))$species ) ) == 1 ){
        print(unique( (loss_tbl%>%filter(representative==rep))$species ) )
        single <- c(single, rep)
    }
}
loss_tbl <- loss_tbl %>% filter(!(species %in% single))

# pdf("strain_loss_new.pdf", width=3, height=3)

this <- loss_tbl %>% filter(!is.na(coverage_loss), complexity!='low')
# randombee_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, box=TRUE, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', clab='Assembler', colorscheme=material, alpha=0.3) + scale_y_reverse()
box_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', flab='Assembler', colorscheme=material) + scale_y_reverse(limits=c(70,0))

this <- loss_tbl %>% filter(!is.na(coverage_loss), complexity!='low', ANI<1)
# randombee_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, box=TRUE, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', clab='Assembler', colorscheme=material, alpha=0.3) + scale_y_reverse()
box_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', flab='Assembler', colorscheme=material) + scale_y_reverse(limits=c(70,0))

this <- loss_tbl %>% filter(!is.na(coverage_loss), complexity!='low', ANI==1)
# randombee_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, box=TRUE, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', clab='Assembler', colorscheme=material, alpha=0.3) + scale_y_reverse()
box_plot(this, this$abundance_bins, this$coverage_loss, this$assembler, xlab='Fold Genome Coverage', ylab='Loss of Genome Fraction [%]', flab='Assembler', colorscheme=material) + scale_y_reverse(limits=c(70,0))
# dev.off()


### BINNERS AND SINGE SAMPLE COMPARISON ###

comp_order <- c('low', 'medium', 'complex')
the_total_tbl$complexity <- factor(the_total_tbl$complexity, levels = comp_order)
the_total_tbl <- the_total_tbl %>% arrange(desc(complexity))

the_total_tbl <- the_total_tbl %>% mutate(abundance_bins=cut(genome_coverage, breaks=c(-Inf, 1, 2.5, 5, 7.5, 10, +Inf), labels=c('< 1x', '1x - 2.5x', '2.5x - 5x', '5x - 7.5x', '7.5x - 10x', '> 10x'))) 

# pdf("binnes_vs_ss_assembly_500bp_new.pdf", width=6, height=4)

this <- the_total_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '500b')

box_plot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Fold Genome Coverage', ylab='Genome fraction [%]', flab='Assembly Strategy', colorscheme=material)
box_plot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembly Strategy', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombee_plot(this, this$assembler, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembly Strategy', colorscheme=material)
this <- the_total_tbl %>% filter(complexity == 'medium', !(assembler %in% c('LSA_33'))) %>% filter(cutoff == '500b')
randombee_plot(this, this$assembler, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembly Strategy', colorscheme=material)
randombee_plot(this, this$assembler, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembly Strategy', colorscheme=material)
randombee_plot(this, this$assembler, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembly Strategy', colorscheme=material) 

ecdf_plot(this, this$genome_fraction, col=this$assembler, colorscheme=material, clab='Assembly Strategy', xlab='Cumulative Density', ylab='Genome Fraction')
ecdf_plot(this, this$NGA50, col=this$assembler, colorscheme=material, clab='Assembly Strategy', xlab='Cumulative Density', ylab='NGA50')
ecdf_plot(this, this$LGA50, col=this$assembler, colorscheme=material, clab='Assembly Strategy', xlab='Cumulative Density', ylab='LGA50')
ecdf_plot(this, this$misassemblies, col=this$assembler, colorscheme=material, clab='Assembly Strategy', xlab='Cumulative Density', ylab='Misassemblies')

# dev.off()

# pdf("binnes_vs_ss_assembly_10kb_new.pdf", width=6, height=4)

this <- the_total_tbl %>% filter(complexity == 'medium') %>% filter(cutoff == '10kb')

box_plot(this, this$abundance_bins, this$genome_fraction, fill=this$assembler, xlab='Fold Genome Coverage', ylab='Genome fraction [%]', flab='Assembler', colorscheme=material)
box_plot(this, this$assembler, this$genome_fraction, fill=this$assembler, xlab='Assembler', ylab='Genome fraction [%]', flab='Complexity', colorscheme=material)
randombee_plot(this, this$assembler, this$genome_fraction, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Genome fraction [%]', clab='Assembler', colorscheme=material)
randombee_plot(this, this$assembler, this$misassemblies, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='Misassemblies', clab='Assembler', colorscheme=material)
randombee_plot(this, this$assembler, this$NGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='NGA50', clab='Assembler', colorscheme=material)
randombee_plot(this, this$complexity, this$LGA50, color=this$assembler, alpha=0.2, box=T, xlab='Complexity', ylab='LGA50', clab='Assembler', colorscheme=material)

# dev.off()


### EVALUATION OF BINNERS/CLUSTERING ###
strains_part <- strains_tbl %>% filter(complexity=='medium')
binners_count_info <- inner_join(strains_part, binners_count_info) %>% inner_join(taxonomy_tbl, by=c('representative'='X3'))
unique(binners_count_info$family)
binners_count_info$family[binners_count_info$family=='NA Bacteroidales fam. [C Bacteroidaceae/Porphyromonadaceae]'] <- 'Bacteroidales'
binners_count_info$family[binners_count_info$family=='NA Bacteria fam. incertae sedis'] <- 'Incertae sedis'
binners_count_info$family <- str_replace(binners_count_info$family, '\\d+ (\\w+)', '\\1')
status_order <- c('representative', 'strain_0', 'strain_1', 'strain_2', 'strain_3', 'strain_4', 'strain_5', 'strain_6')
binners_count_info$status <- factor(binners_count_info$status, levels=status_order)

binners_count_info <- binners_count_info %>% mutate(rep_spec=paste(as.character(status), as.character(representative), sep='_')) 


# pdf("binners_clustering_SAMEA3136631_new.pdf", width=8, height=5)

title <- 'CONCOCT'
this <- binners_count_info %>% filter(assembler=='CONCOCT', sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count_cluster)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE,  colorscheme=many)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)

bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)



title <- 'LSA 55 Clusters'
this <- binners_count_info %>% filter(assembler=='LSA_55', sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count_cluster)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)

bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)


title <- 'LSA 33 Clusters'
this <- binners_count_info %>% filter(assembler=='LSA_33', sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count_cluster)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)

bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)


title <- 'LSA 15 Clusters'
this <- binners_count_info %>% filter(assembler=='LSA_15', sample %in% c('SAMEA3136631'))
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(read_count_cluster)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$species, xlab='Cluster', ylab='Read Abundance', title=title, flegend=FALSE, colorscheme=many)

bar_plot_count(this, this$cluster, this$read_count_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)
bar_plot_count(this, this$cluster, this$read_abundance_cluster, fill=this$family, xlab='Cluster', ylab='Read Abundance', title=title, colorscheme=many, alpha=this$representative)

# dev.off()











# pdf("binnes_enrichment_new.pdf", width=8, height=5)

this <- binners_count_info %>% filter(assembler=='CONCOCT', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'CONCOCT'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=1, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_55', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 55 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=1, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_33', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 33 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=1, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_15', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 15 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=1, colorscheme=many)

### ALPHA ###

this <- binners_count_info %>% filter(assembler=='CONCOCT', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'CONCOCT'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=0.6, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_55', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 55 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=0.6, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_33', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 33 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=0.6, colorscheme=many)

this <- binners_count_info %>% filter(assembler=='LSA_15', rel_abundance>0, max_rel_abundance>0, read_abundance_cluster==max_rel_abundance)
title <- 'LSA 15 Clusters'
loglog_plot(this, this$rel_abundance, this$max_rel_abundance, this$family, xlab='Background relative abundance', ylab='Maximum relative abundance in any cluster', clab='Family', alpha=0.6, colorscheme=many)

# dev.off()



#### CLUSTER ANALYSIS W/ PCOA ####

# pdf("cluster_bars.pdf", width=8, height=5)

gradscale <- color_gradient_scrambled(end='#012824')

this <- binners_count_across %>% filter(assembler=='LSA_55')
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(count)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance', colorscheme=gradscale)

this <- binners_count_across %>% filter(assembler=='LSA_33')
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(count)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance', colorscheme=gradscale)

this <- binners_count_across %>% filter(assembler=='LSA_15')
this <-  this %>% 
    group_by(cluster) %>% 
    mutate(cluster_sum = sum(count)) %>%
    arrange(desc(cluster_sum))
this$cluster <- factor(this$cluster, levels = unique(this$cluster))
bar_plot_count(this, this$cluster, this$abundance, fill=this$species, xlab='Cluster', ylab='Read Abundance', colorscheme=gradscale)

# dev.off()

# pdf("cluster_pcoa.pdf", width=8, height=5)

this <- lsa_55_jsd
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='55 Clusters')
this <- lsa_33_jsd
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='33 Clusters')
this <- lsa_15_jsd
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='15 Clusters')
this <- lsa_55_euclid
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='55 Clusters Euclidean')
this <- lsa_33_euclid
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='33 Clusters Euclidean')
this <- lsa_15_euclid
dot_plot(this, this$Axis.1, this$Axis.2, this$owner, xlab='PC1', ylab='PC2', clab='Main species in cluster', colorscheme=many, size=3, title='15 Clusters Euclidean')
# dev.off()



library(pheatmap)

# pdf("cluster_heatmap_euclid.pdf", width=10, height=10)
pheatmap(lsa_55_euclid_matrix, color=magma(100))
pheatmap(lsa_33_euclid_matrix, color=magma(100))
pheatmap(lsa_15_euclid_matrix, color=magma(100))
# dev.off()

# pdf("cluster_heatmap_jsd.pdf", width=10, height=10)
pheatmap(lsa_55_jsd_matrix, color=colorRampPalette(c('white', '#550527'))(100))
pheatmap(lsa_33_jsd_matrix, color=colorRampPalette(c('white', '#550527'))(100))
pheatmap(lsa_15_jsd_matrix, color=colorRampPalette(c('white', '#550527'))(100))

pheatmap(lsa_55_jsd_matrix, color=colorRampPalette(c('#550527', 'white'))(100))
pheatmap(lsa_33_jsd_matrix, color=colorRampPalette(c('#550527', 'white'))(100))
pheatmap(lsa_15_jsd_matrix, color=colorRampPalette(c('#550527', 'white'))(100))
# dev.off()


lsa_55_jsd_long <- gather_matrix(lsa_55_jsd_matrix, key = c('cluster_1', 'cluster_2'), value='dist')
lsa_15_jsd_long <- gather_matrix(lsa_15_jsd_matrix, key = c('cluster_1', 'cluster_2'), value='dist')


dendro <- as.dendrogram(hclust(d = dist(x = lsa_33_jsd_matrix)))
ggdendrogram(data = dendro, rotate = TRUE)
cluster_order <- rownames(lsa_33_jsd_matrix)[order.dendrogram(dendro)] %>% print()
lsa_33_jsd_long <- gather_matrix(lsa_33_jsd_matrix, key = c('cluster_1', 'cluster_2'), value='dist')
this <- lsa_33_jsd_long %>% group_by(cluster_1, cluster_2) %>% arrange(desc(dist))
this$cluster_1 <- factor(this$cluster_1, levels=cluster_order)
this$cluster_2 <- factor(this$cluster_2, levels=cluster_order)
heatmap_plot(this, this$cluster_1, this$cluster_2, this$dist)



#### POOLED ASSEMBLY ###

method_order <- c('SAMEA_10', 'SAMEA_30', 'lsa_33', 'concoct_66')
pooled_tbl$method <- factor(pooled_tbl$method, levels=method_order)
cutoff_order <- c('500b', '10kb')
pooled_tbl$cutoff <- factor(pooled_tbl$cutoff, levels=cutoff_order)
pooled_tbl <- pooled_tbl %>% mutate(misassembled_fraction=missassembled_contigs_length/total_length)

material <- c('#2598CC', '#05686C', '#7E2382', '#A10702', '#FFB300')

pdf("pooled_assembly.pdf", width=3.5, height=3)

this <- pooled_tbl

that <- this %>% filter(cutoff=='500b')

median((that %>% filter(method=='lsa_33'))$genome_fraction)
box_plot(that, that$method, that$genome_fraction, fill=that$method, xlab='Cutoff', ylab='Genome fraction [%]', flab='Method', colorscheme=material, ylim=c(50, 100))

randombee_plot(that, that$method, that$genome_fraction, color=that$method, xlab='Cutoff', ylab='Genome fraction [%]', clab='Method', colorscheme=material)

that <- this %>% filter(cutoff=='10kb')
box_plot(that, that$method, that$genome_fraction, fill=that$method, xlab='Cutoff', ylab='Genome fraction [%]', flab='Method', colorscheme=material)

this <- pooled_tbl %>% filter(method=='lsa_33')
this$LGA50

box_plot(this, this$cutoff, this$misassembled_fraction, fill=this$method, xlab='Cutoff', ylab='misassembled_fraction', flab='Method', colorscheme=material, ylim=c(0,0.4))

box_plot(this, this$method, this$NGA50, fill=this$method, xlab='Cutoff', ylab='NGA50', flab='Method', colorscheme=material)
box_plot(this, this$method, this$LGA50, fill=this$method, xlab='Cutoff', ylab='LGA50', flab='Method', colorscheme=material, ylim=c(0, 800))

box_plot(this, this$cutoff, this$largest_contig, fill=this$method, xlab='Cutoff', ylab='largest_contig', flab='Method', colorscheme=material)
box_plot(this, this$cutoff, this$total_length, fill=this$method, xlab='Cutoff', ylab='total_length', flab='Method', colorscheme=material, ylim=c(0, 5e+7))

dev.off()

box_plot(this, this$method, this$total_length, fill=this$method, xlab='Cutoff', ylab='total_length', flab='Method', colorscheme=material)
box_plot(this, this$method, this$total_length_1000, fill=this$method, xlab='Cutoff', ylab='total_length_1000', flab='Method', colorscheme=material)
box_plot(this, this$method, this$total_length_10000, fill=this$method, xlab='Cutoff', ylab='total_length_10000', flab='Method', colorscheme=material)
box_plot(this, this$method, this$total_length_50000, fill=this$method, xlab='Cutoff', ylab='total_length_50000', flab='Method', colorscheme=material)




