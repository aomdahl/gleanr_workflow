#!/bin/bash

###March 29, 2024
# Helper script to submit all extractions of GWAS on Rockfish. Yay.
echo -e '#!/bin/bash' > bash_default.tmp
echo "#SBATCH -p express" >> bash_default.tmp
echo "#SBATCH -J partial_extraction" >> bash_default.tmp
echo "#SBATCH --time=150:00" >> bash_default.tmp
echo "source /data/apps/go.sh ### for safety reasons" >> bash_default.tmp
echo "set -e" >> bash_default.tmp
echo "cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts" >> bash_default.tmp

cp bash_default.tmp slurm_first_quarter.sh
sed -i 's/partial_extraction/first_extraction/g' slurm_first_quarter.sh
ls /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/* | head -n 75 | sed 's/^/bash /g' >> slurm_first_quarter.sh

cp bash_default.tmp slurm_second_quarter.sh
sed -i 's/partial_extraction/second_extraction/g' slurm_second_quarter.sh
ls /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/* | head -n 150 | tail -n 75 | sed 's/^/bash /g' >> slurm_second_quarter.sh

cp bash_default.tmp slurm_third_quarter.sh
sed -i 's/partial_extraction/third_extraction/g' slurm_third_quarter.sh
ls /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/* | head -n 230 | tail -n 80 | sed 's/^/bash /g' >> slurm_third_quarter.sh

cp bash_default.tmp slurm_fourth_quarter.sh
sed -i 's/partial_extraction/fourth_extraction/g' slurm_fourth_quarter.sh
ls /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/extract_EUR_scripts/* | tail -n 87  | sed 's/^/bash /g' >> slurm_fourth_quarter.sh

mkdir -p /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files/EUR_only/
sbatch slurm_first_quarter.sh
sbatch slurm_second_quarter.sh
sbatch slurm_third_quarter.sh
sbatch slurm_fourth_quarter.sh
