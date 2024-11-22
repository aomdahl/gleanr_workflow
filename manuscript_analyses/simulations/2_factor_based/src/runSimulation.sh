#!/bin/bash
#To run a simulation, create your desired yaml file and specify an output directory

# Function to display usage
usage() {
  echo "Usage: $0 [-g|-a|-f] <yaml> <output_dir>"
  echo "  -g  Only generate the GWAS data for simulation"
  echo "  -a  Run the entire script"
  echo "  -s  Run the scoring bit at the end" 
  echo "  -f  Run only the matrix factorization step"
  echo "  Note: -g -g and -a may be used in combination if you wish to generate GWAS and factorize, but not evaluate. "
  exit 1
}


# Parse command-line options
RUN_GWAS=false
RUN_ALL=false
RUN_FACTORIZATION=false
RUN_SCORING=false
while getopts "gafs" opt; do
  case $opt in
    g) RUN_GWAS=true ;;
    a) RUN_ALL=true ;;
    f) RUN_FACTORIZATION=true ;;
    s) RUN_SCORING=true ;;
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

if [ "$RUN_GWAS" = false ] && [ "$RUN_FACTORIZATION" = false ] && [ "$RUN_ALL" = false ] &&  [ "$RUN_SCORING" = false ]; then
    echo "Nothing specified to run."
    echo "Program will now conclude"
    exit 0
fi


echo "Currently running $b"

#First: create all the simulations. with niter, each at a different seed
NITER=$(grep "iter" $YAML | cut -f 2 -d ",")
if [ "$RUN_GWAS" = true ] || [ "$RUN_ALL" = true ]; then
  for i in $(eval echo "{1..$NITER}"); do
    echo "iter $i"
    Rscript src/generateGWASFactors.R --input $YAML -o ${ODIR}/sim${i} --seed ${i}
  done
fi

if [ "$RUN_GWAS" = true ] && [ "$RUN_FACTORIZATION" = false ]; then
    echo "GWAS for simulations have been generated."
    echo "Program will now conclude"
    exit 0
fi


# Run matrix factorization
if [ "$RUN_FACTORIZATION" = true ] || [ "$RUN_ALL" = true ]; then
  mkdir -p ${ODIR}/factorization_results
  MNAMES=$(grep "test_methods" $YAML | cut -f 2 -d "," | sed 's/:/,/g')
  K=$(grep -e "^K," $YAML | cut -f 2 -d ",")
  BIC=$(grep "bic_param," $YAML | cut -f 2 -d ",")
  INIT=$(grep "init," $YAML | cut -f 2 -d ",")
  SCALE=$(grep "scale," $YAML | cut -f 2 -d ",")
  SHRINK=$(grep "covar_shrinkage," $YAML | cut -f 2 -d ",")
  SCALE_VAR=""
  SHRINK_VAR=""

  if [ "$SCALE" = "TRUE" ]; then
    SCALE_VAR="--step_scaling"
  fi

  if [ -n "$SHRINK" ]; then
    SHRINK_VAR="--WLgamma $SHRINK"
    echo "shrinkage included"
  fi

  for i in $(eval echo "{1..$NITER}"); do
    echo "Simulation iter $i"
    #Rscript /home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/matrix_factorization.R \
    Rscript ../../src/matrix_factorization.R \
     --se_data ${ODIR}/sim${i}.std_error.txt --beta_data ${ODIR}/sim${i}.effect_sizes.txt --seed ${i} \
      --z_scores ${ODIR}/sim${i}.z.txt -n ${ODIR}/sim${i}.N.txt \
      --outdir ${ODIR}/factorization_results/sim${i}. --only_run $MNAMES --K $K --no_plots --bic_var $BIC --init_mat $INIT \
      --C ${ODIR}/sim${i}.c_matrix.txt $SCALE_VAR $SHRINK_VAR
  done
fi

if [ "$RUN_FACTORIZATION" = true ] && [ "$RUN_ALL" = false ] && [ "$RUN_SCORING" = false ]; then
    echo "Factorization has been succcesfully performed."
    echo "Script will now conclude"
    exit 0
fi

# Evaluate simulations
if [ "$RUN_ALL" = true ] || [ "$RUN_SCORING" = true ]; then
  echo "------------------------------------"
  echo "Beginning simulation scoring"
  Rscript src/evaluateSims.R --output ${ODIR}/factorization_results/summary.noscale --plot --yaml $YAML --sim_path ${ODIR}/factorization_results/
  Rscript src/evaluateSims.R --output ${ODIR}/factorization_results/summary --plot --yaml $YAML --sim_path ${ODIR}/factorization_results/ --scale_data

  echo "END OF FULL SIMULATION"
  echo " "
  echo " "
fi
