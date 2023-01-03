#!/bin/bash
SRC=/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/
ID=$1 #what you want it to be called. Must be an aboslute path, or this breaks on cd
FACTORS=$2 #Path to the factor file
VARIANTS="/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5"
hapmap_list="/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt"
mkdir -p $ID
ml anaconda
conda activate renv
set -e

if [[ "$3" == "NAMES" ]]; then
	A=".names"
else
	A=""
fi
echo "full_hapmap3_snps.z$A.tsv"
if [[ "$4" != "" ]]; then
       VARIANTS=$4
fi
#echo $VARIANTS

#Project
echo "Starting with projection"
echo "Rscript $SRC/projectSumStats.R --output $ID/projected_hapmap3_loadings.txt --sumstats $VARIANTS/full_hapmap3_snps.z.tsv --id_type RSID --factors $FACTORS"
Rscript $SRC/projectSumStats.R --output $ID/projected_hapmap3_loadings.txt --sumstats $VARIANTS/full_hapmap3_snps.z$A.tsv --id_type RSID --factors $FACTORS --proj_method meta

#Make the GWAS files
echo "Making the gwas files..."
mkdir -p $ID/loading_ss_files_meta
Rscript /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/buildSumStats.R --projected_loadings $ID/projected_hapmap3_loadings.txt --samp_file $VARIANTS/full_hapmap3_snps.n$A.tsv --hapmap_list $hapmap_list --output $ID/loading_ss_files_meta/ --factors $FACTORS --samp_counts avg

#LDSC projection
bash src/tissue_testing_only.sh $ID $FACTORS

#rule ldsc_visualize:
bash src/visualize_LDSC_custom.sh $ID/F_ldsc_enrichment_Multi_tissue_chromatin/ 

 


