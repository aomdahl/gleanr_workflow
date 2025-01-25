#!/bin/bash
ml anaconda
conda activate renv
set -e
#cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/simulations/2_factor_based
DEFAULT=/data/abattle4/aomdahl1/reference_data/default_slurm_file.sh
#Build the simulations
Rscript src/buildSimDataClean.R

#Script in previous location used this path, no longer.
#This path specifies both where your directory for the project, which includes your code
#NEW_DIR=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/
NEW_DIR=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/simulations/2_factor_based/
OUTDIR=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/
buildSims()
{
  TYPE=$1
  SEED_APPEND=$2
  echo "$TYPE"
  #Build the appropriate directories
  mkdir -p ${OUTDIR}/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/${TYPE}/
  mkdir -p ${OUTDIR}/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/${TYPE}/

  #Rscript src/BuildYamlFiles.R simulating_factors/custom_easy/setting_files/final_sim_${TYPE}.csv simulating_factors/custom_easy/yaml_files/final_sims_june_2024/${TYPE}/ -c simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/${TYPE}/ > run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt

  Rscript src/BuildYamlFiles.R --param_file ${NEW_DIR}/setting_files/final_sim_${TYPE}.csv --output_path ${OUTDIR}/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/${TYPE}/ \
	  -s ${SEED_APPEND}  --commands ${OUTDIR}simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/${TYPE}/ > ${OUTDIR}/run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt

  #Determine how to split commands by number (help from chatgpt)
  NUM_CALLS=$(wc -l ${NEW_DIR}/setting_files/final_sim_${TYPE}.csv | cut -f 1 -d " ")
  lines_per_file=$(( NUM_CALLS / 3 ))   # Integer division to get base lines per file:
  remainder=$(( NUM_CALLS % 3 ))  # Calculate remainder
  lines_for_file1=$lines_per_file
  inter_file2=$(( NUM_CALLS -lines_for_file1 ))
  lines_for_file2=$lines_per_file
  lines_for_file3=$(( lines_per_file + remainder ))   # For file 3, add the remainder (so we don't lose those lines):

  #Previously hardcoded in numbers, 20,20,15. Now updated, so the current command #'s might not aline
  cat $DEFAULT <(echo "cd ${NEW_DIR}") \
    <(head -n ${lines_for_file1} ${OUTDIR}/run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt) | sed 's/time=5:0:0/time=20:0:0/g' \
    | sed 's/runSimulation.sh/runSimulation.sh -a/g' | sed "s/index_genotype_files/${TYPE}_part1/g" \
    > ${NEW_DIR}/run_scripts/${TYPE}_gwas.first_third.sh
    #> run_scripts/final_sim_june_2024/${TYPE}_gwas.first_third.sh

  #Second set
  cat $DEFAULT <(echo "cd ${NEW_DIR}") \
  <(tail -n ${inter_file2} ${OUTDIR}/run_scripts/final_sim_june_2024/build_${TYPE}_gwas.commands.txt | head -n ${lines_for_file2}) \
    | sed 's/time=5:0:0/time=20:0:0/g' | sed 's/runSimulation.sh/runSimulation.sh -a/g' \
    | sed "s/index_genotype_files/${TYPE}_part2/g" \
    > ${NEW_DIR}/run_scripts/${TYPE}_gwas.second_third.sh
    #> run_scripts/final_sim_june_2024/${TYPE}_gwas.second_third.sh

  #Last one
  cat "$DEFAULT" \
    <(echo "cd ${NEW_DIR}") \
    <(tail -n "${lines_for_file3}" ${OUTDIR}/run_scripts/final_sim_june_2024/build_"${TYPE}"_gwas.commands.txt) \
  | sed 's/time=5:0:0/time=15:0:0/g' \
  | sed 's/runSimulation.sh/runSimulation.sh -a/g' \
  | sed "s/index_genotype_files/${TYPE}_part3/g" \
  > ${NEW_DIR}/run_scripts/"${TYPE}"_gwas.third_third.sh
  #> run_scripts/final_sim_june_2024/"${TYPE}"_gwas.third_third.sh


  echo "All necessary files built. Please inspect"
}

#Main simulations
#buildSims no_overlap 1000
#buildSims 1b_overlap 2000
#buildSims 2b_overlap 3000

#Special cases- no dense factors, muliple dense factors, and some "dens-er" (mid-density) factors

#buildSims special_2b_overlap 4000
#buildSims special_1b_overlap 5000
#buildSims special_no_overlap 6000
buildSims dense_2b_overlap 7000
buildSims dense_1b_overlap 8000
buildSims dense_no_overlap 9000




#cd run_scripts/final_sim_june_2024/
cd ${NEW_DIR}/run_scripts/
for i in *_third.sh; do
	echo $i
	#sbatch $i;
done
