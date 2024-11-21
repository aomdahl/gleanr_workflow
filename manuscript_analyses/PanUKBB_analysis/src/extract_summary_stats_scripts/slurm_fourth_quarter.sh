#!/bin/bash
#SBATCH -p express
#SBATCH -J fourth_extraction
#SBATCH --time=120:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-whr-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-C34-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-C44-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-D12-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-G56-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-I48-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-M16-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-M75-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-N20-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/icd10-N81-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-174-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-250.2-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-300-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-317-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-361-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-365-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-366-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-411.4-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-454.1-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-455-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-471-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-495-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-550.1-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-550.4-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-550-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-562-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-574.1-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-728.71-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-735.3-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-788-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-835-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/phecode-960.2-both_sexes.tsv.extract.sh
