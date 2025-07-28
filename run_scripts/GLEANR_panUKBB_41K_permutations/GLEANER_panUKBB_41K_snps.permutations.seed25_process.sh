#!/bin/bash
#SBATCH --job-name=seed25
#SBATCH --time=3:0:0
#SBATCH --partition=shared
#SBATCH  --mem=10G
source /data/apps/go.sh ### for safety reasons
set -e
ml anaconda
conda activate renv
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization
BD=/scratch16/abattle4/ashton/snp_networks/
mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/assess_test_seed25_100/
Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed25_/permute_100.RData \
--gene_sets ${BD}/custom_l1_factorization/manuscript_analyses/PanUKBB_analysis/analysis/platelet_diseases_gene_lists.txt \
--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/assess_test_seed25_100/ \
--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData -t

#mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/assess_test_seed25_200/
#Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed25_/permute_200.RData \
#--gene_sets ${BD}/custom_l1_factorization/manuscript_analyses/PanUKBB_analysis/analysis/platelet_diseases_gene_lists.txt \
#--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/assess_test_seed25_200/ \
#--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData
