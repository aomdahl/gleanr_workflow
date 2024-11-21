#!/bin/bash
#SBATCH -p shared
#SBATCH -J premunge_all
#SBATCH --time=7:00:00
source /data/apps/go.sh ### for safety reasons
set -e
ml anaconda
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization
  REFSNPS="/data/abattle4/aomdahl1/reference_data/PanUKBB_EUR_LDSC_SNPS.tsv"
  zcat /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_ldsc_files/UKBB.ALL.ldscore/UKBB.EUR.l2.ldscore.gz | awk '(NR > 1){print $1":"$3" "$2}' > $REFSNPS
  #zcat /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_ldsc_files/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz | awk '(NR > 1){print $1":"$3" "$2}' > $REFSNPS
  #Note that these SNPs, as described in https://pan-dev.ukbb.broadinstitute.org/docs/ld, are high quality HM3 variants in autosomes, not in MHC, biallelic, INFOO > 0.9, MAF > 0.01 and gnomad
  #Not sure I totally believe this because I ended up with some rare variants in my downstream analysis based on these SNPs. Hmm. Also that would exclude some.
  trait_file="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/panUKBB_complete.studies.tsv"
  ODIR="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/"
  REFSNPS="/data/abattle4/aomdahl1/reference_data/PanUKBB_EUR_LDSC_SNPS.tsv"
  mkdir -p $ODIR
  conda activate std
  python src/munge_precursor.py --study_list $trait_file  --output $ODIR --mod --keep_maf --merge_alleles "" --rsid_ref $REFSNPS

