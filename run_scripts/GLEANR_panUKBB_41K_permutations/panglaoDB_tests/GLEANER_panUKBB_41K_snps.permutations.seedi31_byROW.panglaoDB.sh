#!/bin/bash
#SBATCH --job-name=assess_test_seed31
#SBATCH --time=3:0:0
#SBATCH --partition=shared
#SBATCH  --mem=10G
source /data/apps/go.sh ### for safety reasons
set -e
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization
ml anaconda
conda activate renv

BD=/scratch16/abattle4/ashton/snp_networks/
PG_GENES=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/panglao_genes.txt
mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed31_100_byROW/
Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed31_byROW/permute_100.RData \
--gene_sets ${PG_GENES} \
--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed31_100_byROW/ \
--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData


mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed31_200_byROW/
Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed31_byROW/permute_200.RData \
--gene_sets ${PG_GENES} \
--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed31_200_byROW/ \
--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData

