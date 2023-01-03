#!/bin/bash
SRC=/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/
ID=$1 #what you want it to be called. Must be an aboslute path, or this breaks on cd
LDSC_REF_DAT="/data/abattle4/aomdahl1/reference_data/ldsc_reference"
echo "Note: default setting is to check if file is already made"
ml anaconda
conda activate py27
pwd
cd $LDSC_REF_DAT
mkdir -p $ID/meta_ldsc_enrichment_Multi_tissue_chromatin/
for l in $ID/loading_ss_files_F/*.gz; do
        
	echo $l 
        T=`basename $l`
        FACTOR=${T%.gz}
        echo $FACTOR
        if test -f "$ID/meta_ldsc_enrichment_Multi_tissue_chromatin/$FACTOR.multi_tissue.cell_type_results.txt"; then
		echo "$FACTOR.multi_tissue.cell_type_results.txt exists, skipping"
		continue
	fi
	python ~/.bin/ldsc/ldsc.py \
        --h2-cts $l \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out $ID/meta_ldsc_enrichment_Multi_tissue_chromatin/$FACTOR.multi_tissue \
        --ref-ld-chr-cts Multi_tissue_chromatin.ldcts \
        --w-ld-chr weights_hm3_no_hla/weights. 
done


