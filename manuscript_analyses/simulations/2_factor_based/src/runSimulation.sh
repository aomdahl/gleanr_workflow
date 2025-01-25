#!/bin/bash
#To run a simulation, create your desired yaml file and specify an output directory

# Function to display usage
usage() {
  echo "Usage: $0 [ -g | -a | -f | -s | -c ] <yaml> <output_dir>"
  echo ""
  echo "Options:"
  echo "  -g  Only generate the GWAS data for simulation"
  echo "  -a  Run the entire script (generate GWAS, factorize, evaluate, etc.)"
  echo "  -f  Run only the matrix factorization step"
  echo "  -s  Run only the scoring step at the end"
  echo "  -c  Run only the consistency scoring step"
  echo ""
  echo "Examples of combining options:"
  echo "  -g -f   Generate GWAS and then do factorization, but skip evaluation."
  echo "  -g -s   Generate GWAS and run scoring if data is already factorized, etc."
  echo ""
  echo "Notes:"
  echo "  If you use -a (run all), there's no need to combine it with other flags."
  exit 1
}

# Parse command-line options
RUN_GWAS=false
RUN_ALL=false
RUN_FACTORIZATION=false
RUN_SCORING=false
RUN_CONSISTENCY=false
while getopts "gafsc" opt; do
  case $opt in
    g) RUN_GWAS=true ;;
    a) RUN_ALL=true ;;
    f) RUN_FACTORIZATION=true ;;
    s) RUN_SCORING=true ;;
    c) RUN_CONSISTENCY=true ;;
    *) usage ;;
  esac
done

# Shift the parsed options away
shift $((OPTIND -1))

# Check for required arguments
if [ $# -ne 2 ]; then
  echo "Incorrectly passed in $# arguments"
  usage
fi

YAML=$1
ODIR=$2

# Activate the conda environment
set -e
ml anaconda
conda activate renv

# Create output directory
mkdir -p $ODIR
b=`basename $YAML`

# If no specific actions were selected, exit
if [ "$RUN_GWAS" = false ] && \
   [ "$RUN_FACTORIZATION" = false ] && \
   [ "$RUN_ALL" = false ] && \
   [ "$RUN_SCORING" = false ] && \
   [ "$RUN_CONSISTENCY" = false ]; then
    echo "No actions specified to run."
    echo "Program will now conclude."
    exit 0
fi

echo "Currently running $b"

###############################
# 1. Create all the simulations. with niter, each at a different seed
NITER=$(grep "iter" $YAML | cut -f 2 -d ",")
if [ "$RUN_GWAS" = true ] || [ "$RUN_ALL" = true ]; then
  for i in $(eval echo "{1..$NITER}"); do
    echo "Iteration $i"
    Rscript src/generateGWASFactors.R --input $YAML -o ${ODIR}/sim${i} --seed ${i}
  done
fi

# If only running simulation generation, exit
if [ "$RUN_GWAS" = true ] && \
   [ "$RUN_FACTORIZATION" = false ] && \
   [ "$RUN_ALL" = false ] && \
   [ "$RUN_SCORING" = false ] && \
   [ "$RUN_CONSISTENCY" = false ]; then
    echo "GWAS for simulations have been generated."
    echo "Program will now conclude."
    exit 0
fi

###########################
# 2. Run matrix factorization
if [ "$RUN_FACTORIZATION" = true ] || [ "$RUN_ALL" = true ]; then
  mkdir -p ${ODIR}/factorization_results
  MNAMES=$(grep "test_methods" $YAML | cut -f 2 -d "," | sed 's/:/,/g')
  K=$(grep -e "^K," $YAML | cut -f 2 -d ",")
  BIC=$(grep "bic_param," $YAML | cut -f 2 -d ",")
  INIT=$(grep "^init," $YAML | cut -f 2 -d ",")
  SCALE=$(grep "scale," $YAML | cut -f 2 -d ",")
  SHRINK=$(grep "covar_shrinkage," $YAML | cut -f 2 -d ",")
  SCALE_VAR=""
  SHRINK_VAR=""

  if [ "$SCALE" = "TRUE" ]; then
    SCALE_VAR="--step_scaling"
  fi

  if [ -n "$SHRINK" ]; then
    SHRINK_VAR="--WLgamma $SHRINK"
    #echo "shrinkage included"
  fi

  for i in $(eval echo "{1..$NITER}"); do
    echo "Simulation iter $i"
    #Rscript /home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/matrix_factorization.R \
    Rscript ../src/matrix_factorization.R \
     --se_data ${ODIR}/sim${i}.std_error.txt --beta_data ${ODIR}/sim${i}.effect_sizes.txt --seed ${i} \
      --z_scores ${ODIR}/sim${i}.z.txt -n ${ODIR}/sim${i}.N.txt \
      --outdir ${ODIR}/factorization_results/sim${i}. --only_run $MNAMES --K $K --no_plots --bic_var $BIC --init_mat $INIT \
      --C ${ODIR}/sim${i}.c_matrix.txt $SCALE_VAR $SHRINK_VAR
  done
fi

# If only factorizing, exit
if [ "$RUN_FACTORIZATION" = true ] && \
   [ "$RUN_ALL" = false ] && \
   [ "$RUN_SCORING" = false ] && \
   [ "$RUN_CONSISTENCY" = false ]; then
    echo "Factorization has been successfully performed."
    echo "Script will now conclude."
    exit 0
fi
###########################################
# 3. Check for factorization results before proceeding
if [ "$RUN_ALL" = true ] || [ "$RUN_SCORING" = true ] || [ "$RUN_CONSISTENCY" = true ]; then
  # Check for factorization outputs before running the next steps.
  if [ -z "$(ls -A "${ODIR}/factorization_results" 2>/dev/null)" ]; then
    echo "Error: No factorization outputs found in ${ODIR}/factorization_results."
    echo "Cannot evaluate simulation performance without factorization results. Exiting..."
    exit 1
  fi

fi

###########################################
# 4. Evaluate simulations
if [ "$RUN_ALL" = true ] || [ "$RUN_SCORING" = true ]; then
  echo "------------------------------------"
  echo "Beginning simulation scoring"
  Rscript src/evaluateSims.R --output ${ODIR}/factorization_results/summary.noscale --plot --yaml $YAML --sim_path ${ODIR}/factorization_results/
  Rscript src/evaluateSims.R --output ${ODIR}/factorization_results/summary --plot --yaml $YAML --sim_path ${ODIR}/factorization_results/ --scale_data

  echo "END OF FULL SIMULATION"
  echo " "
  echo " "
fi

###########################################
# 5. Evaluate consistency across simulations
if [ "$RUN_ALL" = true ] || [ "$RUN_CONSISTENCY" = true ]; then
  echo "------------------------------------"
  echo "Beginning simulation consistency scoring"
  Rscript src/consistency_per_prediction.R \
    --output "${ODIR}/factorization_results/summary" \
    --yaml "$YAML" \
    --sim_path "${ODIR}/factorization_results/" \
    --scale_data

  echo "Completed evaluation of simulation consistency."
  echo ""
  echo ""
fi

