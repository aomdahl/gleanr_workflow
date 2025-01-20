#!/bin/bash
set -e
source /data/apps/go.sh ###
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
ml anaconda; conda activate renv
FP=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
P=${FP}/gwas_extracts/panUKBB_complete/
mkdir -p ${FP}/results/panUKBB_complete_41K_sklearn_eBIC/
#Trait list file:
cut -d " " -f 2 trait_selections/panUKBB_complete.studies.tsv > ${P}/panUKBB_complete.trait_list.tsv
Rscript src/gleaner_grid_slurm.R  --gwas_effects ${P}/panUKBB_complete_clumped_r2-0.2.beta.tsv --uncertainty  ${P}/panUKBB_complete_clumped_r2-0.2.se.tsv --trait_names ${P}/panUKBB_complete.trait_list.tsv \
        --covar_matrix ${FP}/ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv \
        --outdir ${FP}/results/panUKBB_complete_41K_sklearn_eBIC/ --fixed_first --nfactors GRID --sample_sd  ${P}/sample_sd_report.tsv --WLgamma Strimmer \
	--covar_se_matrix ${FP}/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv --bic_var sklearn_eBIC  --converged_obj_change 0.001 \
	--job_name final_gleaner_attempt --time 24:00:00 --memory 75G --task BUILD_JOBS
