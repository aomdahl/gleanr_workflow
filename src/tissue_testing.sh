#!/bin/bash
SRC=/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/
ID=$1 #what you want it to be called. Must be an aboslute path, or this breaks on cd
FACTORS=$2 #Path to the factor file
VARIANTS="/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5"
hapmap_list="/work-zfs/abattle4/ashton/reference_data/hapmap_chr_ids.txt"
mkdir -p $ID
ml gcc/5.5.0
ml R
ml python/3.7-anaconda
set -e
#Project
Rscript $SRC/projectSumStats.R --output $ID/projected_hapmap3_loadings.txt --factors $FACTORS --sumstats $VARIANTS/full_hapmap3_snps.z.tsv --id_type RSID
#
#Make the GWAS files
mkdir -p $ID/loading_ss_files
Rscript $SRC/buildSumStats.R --projected_loadings $ID/projected_hapmap3_loadings.txt --samp_file $VARIANTS/full_hapmap3_snps.n.tsv --hapmap_list $hapmap_list --output $ID/loading_ss_files/ 
#
##Get the LDSC reference information
#mkdir -p ldsc_reference
cd ldsc_reference
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
#tar -xvzf Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
#tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
#tar -xvzf weights_hm3_no_hla.tgz
mkdir -p $ID/ldsc_enrichment
source activate ldsc
pwd
for l in $ID/loading_ss_files/*.gz; do
        echo $l 
        T=`basename $l`
        FACTOR=${T%.gz}
        echo $FACTOR
        python /work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
        --h2-cts $l \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out $ID/ldsc_enrichment_${FACTOR}_ \
        --ref-ld-chr-cts Multi_tissue_gene_expr.ldcts \
        --w-ld-chr weights_hm3_no_hla/weights. 
done


#rule ldsc_visualize:

            Rscript $SRC/visualizeLDSC.R --input_dir $ID/ --plot_type "fdr_sig" --output $ID/fdr_heatmap.png"  --extension "*.cell_type_results.txt
            Rscript $SRC/visualizeLDSC.R --input_dir $ID/ --plot_type "horizontal" --output $ID/full_heatmap.png"  --extension "*.cell_type_results.txt

 


