#!/bin/bash


cd ~/Documents/Msc-Biotechnologie/masterarbeit-zeller

# cat results/quast-metrics.txt | while read METRIC; do Rscript nile/scripts/plotting/aggregate_quast_samples.R `ls results/spades/spades_medium_SAMEA31366*/quast-500b/summary/TSV/${METRIC}.tsv` results/spades/spades-medium-500b/agg_500b_${METRIC}.tsv; done

# Rscript nile/scripts/plotting/aggregate_quast_metrics.R `ls results/spades/spades-complex-10kb/agg_*.tsv` results/spades/spades-complex-10kb/spades_complex_10kb.tsv

Rscript nile/scripts/plotting/aggregate_abundance.R `ls nile/data/specs/abundance-complex/*` abundance_complex.tsv