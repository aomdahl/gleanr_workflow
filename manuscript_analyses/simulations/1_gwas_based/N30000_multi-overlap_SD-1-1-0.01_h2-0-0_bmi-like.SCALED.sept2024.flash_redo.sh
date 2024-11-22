#!/bin/bash
#SBATCH -p shared
#SBATCH -J sims_null_low 
#SBATCH --time=10:00:00
#SBATCH --mem=30G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
set -e
ODIR=N30000_multi-overlap_SD-1_h2-0-0_bmi-like.SCALED_Aug2024/
mkdir -p  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/$ODIR/
echo "Unique case where we modified the terms of phenotype contribution"
echo "manually modified, need to revert to 0.5 and 0.5"
echo "Generating simulations"
echo "------------------------------------------------------------"

echo "Simulation data already made, skipping..."
#Rscript simulation_scripts/gwas_centric_vary_overlap_degree_template.R 23 10 $ODIR  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/ 0 0 1 0,10000,15000,20000,23000,26000,30000 30000
echo "------------------------------------------------------------"
#echo ""
#echo "Simulation data generated. Evaluating now..."

#Originally, ran with kroneker for 25 hours, so slow it only did like 150 runs.
#Just doing flash for speed
#Rscript simulation_scripts/factorize_input_simulations.R --input_dir  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/null_covarfactorization.summaries.flash_kroneker.RData --convergence_criteria 0.005 --run_only "FLASH" 
Rscript simulation_scripts/factorize_input_simulations.R --input_dir  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/${ODIR}/null_covarfactorization.summaries.flash_SE_update.RData --convergence_criteria 0.001 --run_only "FLASH" 


conda deactivate
