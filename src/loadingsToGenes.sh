#!/bin/bash
SU="/work-zfs/abattle4/ashton/snp_utils"
ml gcc/5.5.0
ml R
L=$1 #the loadings file. Should have rownames on it.
N=`showHeader ${1} | wc -l`
SNP_LIST=$2
OD=$3
for i in {1..$N}; do
 Rscript /work-zfs/abattle4/ashton/snp_utils/SNPsToBed.R --snp_list ${SNP_LIST} --loading ${L} --loading_num ${i} --sd 2 --rownames --outdir ${OD}/K${i}.test.bed

  Rscript /work-zfs/abattle4/ashton/snp_utils/SNPsToBed.R --snp_list ${SNP_LIST} --loading ${L} --loading_num ${i} --sd 2 --rownames --outdir ${OD}/K${i}.bg.bed --invert_selection
  bash ${SU}/SNPtoNearestGene.sh ${OD}/K${i}.test.bed  ${OD}/K${i}.test.genes.bed
 bash ${SU}/SNPtoNearestGene.sh ${OD}/K${i}.bg.bed  ${OD}/K${i}.bg.genes.bed
done
