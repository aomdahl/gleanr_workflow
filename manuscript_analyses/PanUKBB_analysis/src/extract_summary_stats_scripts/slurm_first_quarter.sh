#!/bin/bash
#SBATCH -p express
#SBATCH -J first_extraction
#SBATCH --time=120:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30600-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30610-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30620-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30650-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30670-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30680-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30700-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30710-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30720-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30730-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30740-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30750-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30760-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30770-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30810-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30830-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30840-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30850-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30870-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/biomarkers-30890-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-1210-both_sexes-1210.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-1747-both_sexes-1.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20002-both_sexes-1226.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20002-both_sexes-1452.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20002-both_sexes-1466.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20116-both_sexes-0.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20117-both_sexes-1.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-2040-both_sexes-2040.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20502-both_sexes-20502.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-20548-both_sexes-1.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-22126-both_sexes-22126.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-22502-both_sexes-22502.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-2257-both_sexes-2257.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-2463-both_sexes-2463.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/categorical-2473-both_sexes-2473.tsv.extract.sh
