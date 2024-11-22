#!/bin/bash
ml anaconda
conda activate renv
cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration
DEFAULT=/data/abattle4/aomdahl1/reference_data/default_slurm_file.sh
#Build the simulations
Rscript src/buildSimDataClean.R



buildSims()
{
TYPE=$1
echo "$TYPE"
    #Build the appropriate directories
    mkdir -p simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/${TYPE}/
    mkdir -p simulating_factors/custom_easy/yaml_files/final_sims_june_2024/${TYPE}/

Rscript src/BuildYamlFiles.R simulating_factors/custom_easy/setting_files/final_sim_${TYPE}.csv simulating_factors/custom_easy/yaml_files/final_sims_june_2024/${TYPE}/ -c simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/${TYPE}/ > run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt

#First 20 commands
cat $DEFAULT <(echo "cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration") <(head -n 20 run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt) | sed 's/time=5:0:0/time=20:0:0/g' | sed 's/runSimulation.sh/runSimulation.sh -a/g' | sed "s/index_genotype_files/${TYPE}_part1/g" > run_scripts/final_sim_june_2024/${TYPE}_gwas.first_third.sh

#Second 20 commands
cat $DEFAULT <(echo "cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration") <(tail -n 35 run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt | head -n 20) | sed 's/time=5:0:0/time=20:0:0/g' | sed 's/runSimulation.sh/runSimulation.sh -a/g' | sed "s/index_genotype_files/${TYPE}_part2/g" > run_scripts/final_sim_june_2024/${TYPE}_gwas.second_third.sh

#Last 15
cat $DEFAULT <(echo "cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration") <(tail -n 15 run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt) | sed 's/time=5:0:0/time=15:0:0/g' | sed 's/runSimulation.sh/runSimulation.sh -a/g' | sed "s/index_genotype_files/${TYPE}_part3/g" > run_scripts/final_sim_june_2024/${TYPE}_gwas.third_third.sh

echo "All necessary files built. Please inspect"
}

buildSims no_overlap
buildSims 1b_overlap
buildSims 2b_overlap
cd run_scripts/final_sim_june_2024/
for i in *_third.sh; do 
	echo $i
	#sbatch $i;
done
