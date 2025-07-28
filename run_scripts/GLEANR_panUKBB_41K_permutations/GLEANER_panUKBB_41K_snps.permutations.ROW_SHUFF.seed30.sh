#!/bin/bash
#SBATCH -J GLEANER_full_permute_30
#SBATCH --time=3:00:00
#SBATCH --partition=shared
#SBATCH --mem=30G
set -e
source /data/apps/go.sh ###
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
ml anaconda; conda activate renv
FP=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
P=${FP}/gwas_extracts/panUKBB_complete/
Rscript src/gleanr_permute.R  --gwas_effects ${P}/panUKBB_complete_clumped_r2-0.2.beta.tsv --uncertainty  ${P}/panUKBB_complete_clumped_r2-0.2.se.tsv --trait_names ${P}/panUKBB_complete.trait_list.tsv \
        --covar_matrix ${FP}/ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv \
        --outdir /scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed30_byROW --fixed_first --nfactors 136 --sample_sd  ${P}/sample_sd_report.tsv --WLgamma Strimmer \
        --covar_se_matrix ${FP}/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv --bic_var sklearn  --converged_obj_change 0.001 \
        --run_data ${FP}/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData \
        --seed 30 --n_permutations 200 --shuff_type rows
