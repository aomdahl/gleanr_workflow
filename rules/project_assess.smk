# load modules
shell.prefix("ml anaconda; conda activate renv;") 
# configurations
#configfile: "config/config.yaml"
src_path="/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/"
LDSC_REF_DAT="/data/abattle4/aomdahl1/reference_data/ldsc_reference"

#For later:
output_type =["loadings", "factors", "weights"]
#make things better with the checkpoint tutorial:https://evodify.com/snakemake-checkpoint-tutorial/
rule all:
    input:
	"results//infertility/p0.001_FULL/flash_backfit_zscores//ldsc_enrichment_Multi_tissue_chromatin/fdr_0.05_heatmap.png"	
        #"results/infertility/p0.001_FULL/flash_backfit_zscores/projectedF_hapmap3_loadings.txt"


rule hapmap_reference: #get the list of hapmap snps for extraction, omitting HLA region
    input:
    output:
        #"data/hm3_no_hla.txt" #being lay, putting in my version for now...
        "/data/abattle4/aomdahl1/reference_data/hm3_nohla.snpids"
    shell:
        """
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz -p data/
        tar -xvzf data/weights_hm3_no_hla.tgz
        for i in "factorization_data/{identifier}.factors.txt"{{1..22}}; do zcat data/weights_hm3_no_hla/weights.${{i}}.l2.ldscore.gz | tail -n +2 | awk '{{print $1":"$3"\t"$2}}' >> {output}; done
        """

rule project_F:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="results/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projected_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' ")

rule project_LM:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="results/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedLM_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' --proj_method std_lm")

rule build_LD_space:
	input:
	output:
		"test.txt"
	run:
		"""
		ml plink
		#make the following its own rule?
		awk '(FNR == NR) {arr[$2]=$1;next} ($1 in arr) {print arr[$1]}' /data/abattle4/aomdahl1/reference_data/hm3_nohla.snpids allele_merge_list.txt | sort -h > 1000G.snp.ids
		plink --bfile /scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/ref --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --ld-snp-list 1000G.snp.ids   --out {output}
		"""
rule project_LD:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="results/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedLM_hapmap3_loadings.txt"
    params:
    run:
        shell("ml anaconda; conda activate std; python --ld_ref")
#option to make a rule project_Surya:
checkpoint prep_enrichment: #format the outputed factors for enrichment analysis
    input:
        hapmap_list="/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt",
        projections="results/{identifier}/projected_hapmap3_loadings.txt",
        sample_counts="gwas_extracts/{identifier}/full_hapmap3_snps.N.tsv"

    output:
        directory("results/{identifier}/loading_ss_files")
    params:
        "results/{identifier}/loading_ss_files/"
    shell:
        """
            mkdir -p {output}
            Rscript {src_path}/buildSumStats.R --projected_loadings {input.projections} --samp_file {input.sample_counts} --hapmap_list {input.hapmap_list} --output {params} --normal_transform
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
        ss="results/{identifier}/loading_ss_files/{factor}.sumstats.gz",
	tiss_dir=LDSC_REF_DAT+"/{tis_ref}.ldcts"
    output:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.cell_type_results.txt",
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.log"
    params:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue",
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
        --w-ld-chr weights_hm3_no_hla/weights.
        cd ../
        """

def aggregate_factors(wildcards):
    checkpoint_output = checkpoints.prep_enrichment.get(**wildcards).output[0]
    factor_numbers = glob_wildcards(f"{checkpoint_output}/F{{factor}}.sumstats.gz").factor
    print(factor_numbers)
    ldsc_files=expand("results/{{identifier}}/ldsc_enrichment_{{tis_ref}}/F{fn}.multi_tissue.cell_type_results.txt", fn = factor_numbers)
    print(ldsc_files)
    return ldsc_files



rule ldsc_visualize:
    input:
       aggregate_factors
    output:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/full_heatmap.png",
	"results/{identifier}/ldsc_enrichment_{tis_ref}/fdr_0.05_heatmap.png", 
	"results/{identifier}/ldsc_enrichment_{tis_ref}/fdr_0.01_heatmap.png",
	"results/{identifier}/ldsc_enrichment_{tis_ref}/wrapped_z.heatmap.png"
    params:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/"
    shell:
        """
            echo {input}
            Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "facet_wrap" --output {output[3]}
	    Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "fdr_sig" --output {output[1]} --fdr 0.05
            Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "fdr_sig" --output {output[2]} --fdr 0.01
            Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "horizontal" --output {output[0]}
        """

rule factors_assessment:
#This isn't perfect. For a cleaner run of this, try:
# bash src/runOnCustomOnes.sh ./factorization_run_lists/7_k_runlist.txt
#where 7_k_runlist.txt is a list of all of the factorizations to analyze.
#In the future, I would like to have this step nicely snakemaked....
    input: #a bit hacky at the moment, but whatever...
        tiss_dir = "results/{identifier}/ldsc_enrichment_{tis_ref}/",
        trait_names = "/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv",
        trait_ids = "/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", 
        factors= "factorization_data/{identifier}.factors.txt"
    output:  "results/{identifier}/factor_simple_scores.txt"
    params: "results/{identifier}"
    shell:
        """
            echo "Assuming using the seed 2 run...."
            Rscript /scratch16/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/src/factAssessment.R --factors {input.factors} \
                --output {params[0]} --simple --ldsc_reference  ldsc_results/seed2_thres0.9_h2-0.1/ \
                --ldsc_dir {input.tiss_dir} --trait.ids {input.trait_ids} --trait.names {input.trait_names}
        """
