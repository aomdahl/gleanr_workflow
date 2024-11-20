#include: "ldsr_pairwise.smk"
# load modules
shell.prefix("ml anaconda; conda activate renv;") 
# configurations
#configfile: "config/config.yaml"
src_path="/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/"
LDSC_REF_DAT="/data/abattle4/aomdahl1/reference_data/ldsc_reference"
GENE_LIBS=["WikiPathways_2023_Human", 
"Chromosome_Location",
"DisGeNET", 
"GO_Biological_Process_2023", 
"GO_Cellular_Component_2023", 
"GO_Molecular_Function_2023", 
"GWAS_Catalog_2023", 
"KEGG_2021_Human", 
"MSigDB_Hallmark_2020", 
"OMIM_Expanded",
"Azimuth_Cell_Types_2021"]
SNP_GENE_MAP="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/41K_openTargets.csv"
analysis_scripts="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/"
nfactors=58
rule all:
	input: 
		expand("results/panUKBB_complete_61K/META_ldsc_enrichment_Multi_tissue_chromatin/F{num}.multi_tissue.log", num=list(range(1,52))),
		
		#expand("results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/F{num}.multi_tissue.log", num=list(range(1,52))),
		#"results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/factor_global_fdr.heatmap.png"

rule selective_pressure:
	input: "results/{identifier}/_final_dat.RData"
	output: "results/{identifier}/selective_pressure/allFigs.RData"
	params:
		odir="results/{identifier}/selective_pressure/"
	shell:
		"""
			Rscript src/evaluateSelectivePressure.R --factorization {input} --output {params.odir}
		"""	

rule top_snps_by_factor:
    input: "results/{identifier}/_final_dat.RData"
    output: "results/{identifier}/top_elements_by_factor/"

rule top_genes_by_factor:
    input:
        snp_gene_map=SNP_GENE_MAP,
        snp_scores="results/{identifier}/latent.loadings.txt"
    #How did I map these previously? Don't need to do it again....
    output:
        background_set=expand("results/{{identifier}}/top_elements_by_factor/F{F}_background_genes.csv", F=range(1,nfactors+1)),
        fg_set = expand("results/{{identifier}}/top_elements_by_factor/F{F}_genes_in_factor.csv", F=range(1,nfactors+1))
    params:
        outdir="results/{identifier}//top_elements_by_factor/"
    shell:
        """
            conda activate renv
            Rscript {analysis_scripts}/getTopGenes.R --snp_gene_map {input.snp_gene_map} --snp_scores {input.snp_scores} --output {params.outdir}
        """

rule gene_set_enrichment:
    input: 
        gene_set="results/{identifier}/top_elements_by_factor/genes/F{factor}.gene_set.txt",
        background_set="results/{identifier}/top_elements_by_factor/genes/F{factor}.background_set.txt",
    output: 
       "results/{identifier}/top_elements_by_factor/gene_factor_enrichments/F{factor}.{gene_lib}.txt",
    shell:
        """
        conda activate std
        python {analysis_scripts}/gsea_enrichr.py -t {input.gene_set} -b {input.background_set} -o {output} -g {gene_lib}
        """


rule join_enrichment_dat:
    input:
        expand("results/{{identifier}}/top_elements_by_factor/gene_factor_enrichments/F{{factor}}.{gene_list}.txt", gene_list=GENE_LIBS)
    output:
        "results/{identifier}/top_elements_by_factor/gene_factor_enrichments/F{factor}.joined_enrichments.txt"
    shell:
        """
        touch {output}
        for i in {input}; do
            tail -n +2 $i >> {output}
        done
        """

