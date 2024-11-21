#!/bin/bash
#SBATCH -p express
#SBATCH -J third_extraction
#SBATCH --time=120:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30080-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30100-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30110-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30120-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30130-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30150-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30180-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30190-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30300-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-30530-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-3062-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-3064-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-3143-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-3148-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-400-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-404-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-4195-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-4282-both_sexes.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-46-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5085-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-50-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5104-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5110-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5116-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5134-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5157-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5257-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-5262-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-6033-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-DBP-both_sexes-combined_medadj_irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-FEV1FVC-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-LDLC-both_sexes-medadj_irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-MCP-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-NAP-both_sexes-irnt.tsv.extract.sh
bash /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/continuous-PP-both_sexes-combined_medadj_irnt.tsv.extract.sh
