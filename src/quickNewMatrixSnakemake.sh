#!/bin/bash
#This script is for quickly taking a factorization and making it 
#inputs are factor #2 id and #1 path and #3 number of cores for the snakemake
set -e
NCORE=$3
ID=$2
        P=$1
        cp $P ./factorization_data/${ID}.factors.txt
        ml python/3.7-anaconda
        source activate mysnakemake
        mkdir -p gwas_extracts/${ID}/
        cp gwas_extracts/PMA_91_2.3_8-11/* gwas_extracts/${ID}/
        snakemake results/${ID}/ldsc_enrichment_Multi_tissue_chromatin/full_heatmap.png -j $NCORE
