#!/bin/bash
#SBATCH --job-name=seed23_panglao
#SBATCH --time=3:0:0
#SBATCH --partition=shared
#SBATCH  --mem=10G
source /data/apps/go.sh ### for safety reasons
set -e
ml anaconda
conda activate renv
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization

#Rscript /scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/panglao_make_gene_list.R
PG_GENES=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/panglao_genes.txt
BD=/scratch16/abattle4/ashton/snp_networks/
mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed23_100/
Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed23_/permute_100.RData \
--gene_sets ${PG_GENES} \
--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed23_100/ \
--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData

mkdir -p ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed23_200/
Rscript src/evaluate_gleanr_permutations.R --input=${BD}/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed23_/permute_200.RData \
--gene_sets ${PG_GENES} \
--num_columns=4,11,23 --output ${BD}/scratch/manuscript_reviews/permute_testing/panglaoDB/assess_test_seed23_200/ \
--ref_analysis ${BD}/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData
