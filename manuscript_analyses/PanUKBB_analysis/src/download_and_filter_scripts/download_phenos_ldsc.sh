#!/bin/bash
#SBATCH -p express
#SBATCH --time=30:00
set -e
source /data/apps/go.sh
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_ldsc_files/
wget https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz
tar -xvzf UKBB.ALL.ldscore.tar.gz
