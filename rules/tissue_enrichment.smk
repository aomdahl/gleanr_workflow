

shell.prefix("ml anaconda; conda activate renv;") 
# configurations
#configfile: "config/config.yaml"
src_path="/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/"
LDSC_REF_DAT="/data/abattle4/aomdahl1/reference_data/ldsc_reference"
snp_reflist="/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz"
#        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/factor_global_fdr.heatmap.png"
rule all:
    input: expand("results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_gene_expr/F{fn}.multi_tissue.cell_type_results.txt", fn=range(1,59))
	#"results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin//factor_global_fdr.heatmap.png"

checkpoint project_none:
    input:
        hapmap_list="/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt",
        loadings="results/{identifier}/latent.loadings.txt",
        factors="results/{identifier}/latent.factors.txt",
        sample_counts="gwas_extracts/{identifier}/{identifier}.n.tsv", #This isn't right for urrent setup, but close enough,
        sample_se="gwas_extracts/{identifier}/{identifier}.se.tsv"
    output:
        directory("results/{identifier}/loading_ss_files_NONE/")
    shell:
        """
        mkdir -p {output}
        conda activate renv
        echo "running the thing now..."
        Rscript src/buildSumStatsDirect.R --loadings {input.loadings} \
         --factors {input.factors} \
         --output {output} \
         --samp_file {input.sample_counts} \
         --samp_se  {input.sample_se} \
         --snp_reflist {snp_reflist}
         echo "done"
        """ 


rule download_enrichment_refs:
    input:
    output:
        "{LDSC_REF_DAT}/Multi_tissue_gene_expr.ldcts"
    params:
        ldscpath={LDSC_REF_DAT}
    shell:
       	"""
            mkdir -p {params.ldscpath}
            cd {params.ldscpath}
            wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
            wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
            wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
            tar -xvzf Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
            tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
            tar -xvzf weights_hm3_no_hla.tgz
        """

#tis_ref may be either "Multi_tissue_chromatin" or Multi_tissue_gene_expr"
rule ldsc_enrichment: #just run for one, then call on the input.
    input:
        ss="results/{identifier}/loading_ss_files_{projection_style}/{factor}.sumstats.gz",
	    tiss_dir=LDSC_REF_DAT + "/{tis_ref}.ldcts"
    output:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.cell_type_results.txt",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.log"
    params:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.multi_tissue",
	LDSC_REF_DAT
    shell:
        """
	ml anaconda
	conda activate py27
	CWD=`pwd`
        cd {params[1]}
        python /scratch16/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
        --h2-cts $CWD/{input.ss} \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out $CWD/{params[0]} \
        --ref-ld-chr-cts {input.tiss_dir} \
        --w-ld-chr weights_hm3_no_hla/weights. --n-blocks 5000
        cd ../
        """

def aggregate_factors(wildcards):
    checkpoint_output = checkpoints.project_none.get(**wildcards).output[0]
    factor_numbers = glob_wildcards(f"{checkpoint_output}/F{{factor}}.sumstats.gz").factor
    print(factor_numbers)
    ldsc_files=expand("results/{{identifier}}/{{projection_style}}_ldsc_enrichment_{{tis_ref}}/F{fn}.multi_tissue.cell_type_results.txt", fn = factor_numbers)
    print(ldsc_files)
    return ldsc_files



rule ldsc_visualize:
    input:
       aggregate_factors
    output:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/full_heatmap.png",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/fdr_0.05_heatmap.png", 
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/fdr_0.01_heatmap.png",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/wrapped_z.heatmap.png",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/factor_tissue_fdr.heatmap.png",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/factor_global_fdr.heatmap.png"
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/plotting_table.csv"
    params:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/"
    shell:
        """
            echo {input}
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "facet_wrap" --output {output[3]}
	    Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "fdr_sig" --output {output[1]} --fdr 0.05
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "fdr_sig" --output {output[2]} --fdr 0.01
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "horizontal" --output {output[0]}
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "factor_tissue_FDR" --output {output[4]}
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "global_tissue_FDR" --output {output[5]}
            Rscript src/visualizeLDSC.R --input_dir {params} --plot_type "data_tab" --output {output[6]}
        """
