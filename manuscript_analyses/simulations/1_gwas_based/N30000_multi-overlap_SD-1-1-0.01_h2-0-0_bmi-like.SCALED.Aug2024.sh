#!/bin/bash
#SBATCH -p shared
#SBATCH -J sims_null_low 
#SBATCH --time=35:00:00
#SBATCH --mem=30G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/simulations/1_gwas_based
ml anaconda
conda activate renv
set -e
ODIR=N30000_multi-overlap_SD-1_h2-0-0_bmi-like.SCALED_Aug2024/
mkdir -p  results/$ODIR/
echo "Generating simulations"
echo "------------------------------------------------------------"

#echo "Simulation data already made, skipping..."
Rscript src/gwas_centric_vary_overlap_degree_template.R 23 10 $ODIR  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/ 0 0 1 0,10000,15000,20000,23000,26000,30000 30000
echo "------------------------------------------------------------"
echo ""
echo "Simulation data generated. Evaluating now..."

Rscript src/factorize_input_simulations.R --input_dir  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/null_covarfactorization.summaries.RData  --with_covar --convergence_criteria 0.001 


conda deactivate
