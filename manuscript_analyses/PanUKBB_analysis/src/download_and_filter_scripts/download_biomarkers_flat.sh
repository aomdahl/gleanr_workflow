c.....bin.bash.....SBATCH..p.express.....SBATCH...time.60.00...
#!/bin/bash
#SBATCH -p express
#SBATCH --time=60:00
set -e
source /data/apps/go.sh
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30610-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30620-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30630-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30640-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30650-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30670-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30680-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30690-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30700-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30710-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30720-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30730-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30740-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30750-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30760-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30770-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30780-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30810-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30830-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30840-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30850-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30860-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30870-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30880-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30890-both_sexes-irnt.tsv.bgz
