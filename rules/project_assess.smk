#include: "ldsr_pairwise.smk"
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
		expand("results/panUKBB_complete_61K/META_ldsc_enrichment_Multi_tissue_chromatin/F{num}.multi_tissue.log", num=list(range(1,52))),
		
		#expand("results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/F{num}.multi_tissue.log", num=list(range(1,52))),
		#"results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/factor_global_fdr.heatmap.png"

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

rule prepProjectionSpace:
        input:
                "results/{identifier}/latent.factors.txt"
                #add this in later- the directory of z-scoresdirector(loadings, the directory of summary stats
        output:
            "results/{identifier}/munged_paths_traits.txt"
        params:
            "results/{identifier}/",
	    "gwas_extracts/{identifier}/"
        shell:
            """
   	#line bleow should be input 1
            readlink -f {params[1]}/*.sumstats.gz | awk -F "." '{{print $0 "\t" $(NF-2)}}' > {params[0]}/ldsc_extract_order.tmp
            tail -n +2 {params[0]}/latent.factors.txt | cut -f 1 -d " " > {params[0]}/trait.order.tmp

            awk '(FNR == NR) {{arr[$2]=$1;next}} ($1 in arr) {{print arr[$1]"\t"$1}}' {params[0]}/ldsc_extract_order.tmp {params[0]}/trait.order.tmp > {output}

            rm {params[0]}/*.tmp
            """
#This should be replaced by a newer version of this rule
rule projectionSpaceOLD:
        input:
                "results/{identifier}/munged_paths_traits.txt"
        output:
            "gwas_extracts/{identifier}/full_hapmap3_snps.Z.tsv", "gwas_extracts/{identifier}/full_hapmap3_snps.N.tsv"
        params:
            "gwas_extracts/{identifier}/"
        shell:
            """
            ml anaconda
            conda activate renv

            Rscript src/ldscToMatrix.R --filepathlist {input} --outdest {params}/full_hapmap3_snps --feature_list "Z,N"
            """

rule projectionSpace:
        input:
             "gwas_extracts/{identifier}/missingness_report.tsv"
        output:
            "gwas_extracts/{identifier}/full_hapmap3_snps.z.tsv", "gwas_extracts/{identifier}/full_hapmap3_snps.n.tsv"
        params:
            "gwas_extracts/{identifier}/",
	    "gwas_extracts/{identifier}/full_hapmap3_snps"
        shell:
            """
            ml anaconda
            conda activate std

            python src/quickGWASIter.py  --type ldsc_custom --output {params[1]} --gwas_list {input} --snp_list "/data/abattle4/aomdahl1/reference_data/ldsc_tutorial_hm3.no_hla.snplist"  --extension .sumstats.gz --gwas_dir  {params[0]}
            """
rule project_F:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedF_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' ")

rule project_LMscaled:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedLMscaled_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --standardize --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' --proj_method std_lm")
#TODO add the adjustment procedure
rule project_META:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedMETA_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' --proj_method meta")

rule project_LM:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedLM_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' --proj_method std_lm")



rule project_LMwht:
    input:
        #factors="factorization_data/{identifier}.factors.txt", #MAJOR CHANGE HERE. MOVING to be in results category
        factors="results/{identifier}/latent.factors.txt",
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.Z.tsv", #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
	mat="ldsr_results/{identifier}/summary_data/gcov_int.tab.csv"
    output:
        "results/{identifier}/projectedLMwht_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID' --proj_method std_lm --decorrelate {input.mat}")

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
        variants="gwas_extracts/{identifier}/full_hapmap3_snps.Z.tsv" #make a sym link if its in a common dir, so we don't have like a jillion copies of larege files.
    output:
        "results/{identifier}/projectedLD_hapmap3_loadings.txt", "results/{identifier}/done.in.ld.txt"
    params:
    run:
        shell("ml anaconda; conda activate std; python --ld_ref")
#option to make a rule project_Surya:


checkpoint prep_enrichment: #format the outputed factors for enrichment analysis
    input:
        hapmap_list="/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt",
        projections="results/{identifier}/projected{projection_style}_hapmap3_loadings.txt",
        sample_counts="gwas_extracts/{identifier}/full_hapmap3_snps.n.tsv",
	factors="results/{identifier}/latent.factors.txt"
    output:
        directory("results/{identifier}/loading_ss_files_{projection_style}/")
    params:
        "results/{identifier}/loading_ss_files_{projection_style}/"
    shell:
        """
	    echo "assuming the 1st columns of input data are SNP IDs"
            mkdir -p {output}
            Rscript ./src/buildSumStats.R --projected_loadings {input.projections} --samp_file {input.sample_counts} --hapmap_list {input.hapmap_list} --output {params} --normal_transform --factors {input.factors}
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
        --w-ld-chr weights_hm3_no_hla/weights. #--n-blocks 5000
        cd ../
        """
rule ldsc_enrichment_panUKBB_jacknife: #just run for one, then call on the input.
    input:
        ss="results/{identifier}/loading_ss_files_{projection_style}/{factor}.sumstats.gz",
	tiss_dir=LDSC_REF_DAT + "/{tis_ref}.ldcts"
    output:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.PanUKBB-ref.multi_tissue.cell_type_results.txt",
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.PanUKBB-ref.multi_tissue.log"
    params:
        "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/{factor}.PanUKBB-ref.multi_tissue",
	LDSC_REF_DAT
    shell:
        """
	ml anaconda
	conda activate py27
	CWD=`pwd`
        cd {params[1]}
        echo "We are in directory {params[1]}"
	pwd
	python /scratch16/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
        --h2-cts $CWD/{input.ss} \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out $CWD/{params[0]} \
        --ref-ld-chr-cts {input.tiss_dir} \
        --w-ld-chr UKBB.ALL.ldscore/UKBB.EUR \
	--n-blocks 5000
        cd ../
        """

def aggregate_factors(wildcards):
    checkpoint_output = checkpoints.prep_enrichment.get(**wildcards).output[0]
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
        """

rule factors_assessment:
#This isn't perfect. For a cleaner run of this, try:
# bash src/runOnCustomOnes.sh ./factorization_run_lists/7_k_runlist.txt
#where 7_k_runlist.txt is a list of all of the factorizations to analyze.
#In the future, I would like to have this step nicely snakemaked....
    input: #a bit hacky at the moment, but whatever...
        tiss_dir = "results/{identifier}/{projection_style}_ldsc_enrichment_{tis_ref}/",
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
