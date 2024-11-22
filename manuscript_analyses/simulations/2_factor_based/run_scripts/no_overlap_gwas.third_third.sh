#!/bin/bash
#SBATCH --job-name=no_overlap_part3
#SBATCH --time=15:0:0
#SBATCH --partition=shared
#SBATCH  --mem=10G
source /data/apps/go.sh ### for safety reasons 
set -e
cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V101_U102_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U102_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V102_U102_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V102_U102_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V103_U102_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V103_U102_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V101_U103_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U103_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V102_U103_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V102_U103_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V103_U103_MAF-mix_eur_N-10000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V103_U103_MAF-mix_eur_N-10000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V101_U101_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U101_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V102_U101_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V102_U101_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V103_U101_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V103_U101_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V101_U102_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U102_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V102_U102_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V102_U102_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V103_U102_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V103_U102_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V101_U103_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U103_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V102_U103_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V102_U103_MAF-mix_eur_N-5000_RHO-none_No-none/ 
bash src/runSimulation.sh -f -s  simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap//V103_U103_MAF-mix_eur_N-5000_RHO-none_No-none.yml simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V103_U103_MAF-mix_eur_N-5000_RHO-none_No-none/ 
