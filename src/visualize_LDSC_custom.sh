#!/bin/bash
SRC=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/
ID=$1 #what you want it to be called. Must be an aboslute path, or this breaks on cd
ml anaconda
conda activate renv
echo $ID
            Rscript $SRC/visualizeLDSC.R --input_dir $ID --plot_type "factor_tissue_FDR" --output $ID/factor_tissue_fdr.heatmap.png  --extension "*.cell_type_results.txt" --fdr 0.05
            Rscript $SRC/visualizeLDSC.R --input_dir $ID --plot_type "facet_wrap" --output $ID/wrapped_z.heatmap.png  --extension "*.cell_type_results.txt"

 


