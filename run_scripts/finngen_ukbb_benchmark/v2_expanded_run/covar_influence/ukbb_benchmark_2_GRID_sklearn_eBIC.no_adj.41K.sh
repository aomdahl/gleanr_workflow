#!/bin/bash
#SBATCH -p shared
#SBATCH -J sklearn_eBIC_41
#SBATCH --time=30:00
#SBATCH --mem=10G
#SBATCH --output=R-%x.%j.out
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
ml anaconda; conda activate renv
Rscript src/gleaner_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2_conservative.beta.tsv --uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//ukbb_benchmark_2_conservative.se.tsv --trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//pheno.names.txt --outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence//BIC-sklearn_eBIC_K-41_no_covar --fixed_first -K "41" --bic_var "sklearn_eBIC" -v 1 \
