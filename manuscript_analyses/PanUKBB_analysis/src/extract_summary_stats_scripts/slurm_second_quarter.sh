#!/bin/bash
#SBATCH -p express
#SBATCH -J second_extraction
#SBATCH --time=120:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-2644-both_sexes-2644.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-3606-both_sexes-3606.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-3799-both_sexes-3799.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-4598-both_sexes-4598.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6148-both_sexes-100.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6149-both_sexes-1.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6149-both_sexes-3.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6149-both_sexes-4.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6149-both_sexes-6.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-6159-both_sexes-1.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-102-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1160-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1180-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1200-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1220-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1687-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-1717-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20016-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20022-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20127-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20154-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20409-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20414-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-20460-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-2139-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-2149-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-22677-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-23100-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-23101-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-23106-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30000-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30010-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30040-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30060-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30070-both_sexes-irnt.tsv.extract.sh
