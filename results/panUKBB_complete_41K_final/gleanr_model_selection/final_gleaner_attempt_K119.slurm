#!/bin/bash
#SBATCH --job-name=final_gleaner_attempt_K119
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --mem=75G
#SBATCH --output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC//%j.%x.out
#SBATCH --error=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC//%j.%x.err
set -e
mkdir -p /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC/
source /data/apps/go.sh 
 cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
ml anaconda
conda activate renv
Rscript src/gleaner_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete_clumped_r2-0.2.beta.tsv --uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete_clumped_r2-0.2.se.tsv --trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv --weighting_scheme B_SE --covar_matrix /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv --converged_obj_change 0.001 --niter 300 --outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K119 --fixed_first --nfactors 119 --bic_var sklearn_eBIC --param_conv_criteria BIC.change --verbosity 0 --ncores 1 --sample_sd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//sample_sd_report.tsv --block_covar 0.2 --covar_se_matrix /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv --subset_seed -1 --WLgamma Strimmer --model_select

