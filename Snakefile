# load modules
shell.prefix("module load gcc/5.5.0; module load R;") 
# configurations
configfile: "config/config.yaml"
src_path="/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/"
#For later:
output_type =["loadings", "factors", "weights"]
factorization_type=["flashr", "pca"]
#make things better with the checkpoint tutorial:https://evodify.com/snakemake-checkpoint-tutorial/
rule all:
    input:
       # expand("results/seed2_thresh0.9_h2-0.1_vars0.01/ldsc_enrichment/F{factor}.multi_tissue.cell_type_results.txt", factor=["1","2","3","4","5","6","7", "8", "9", "10", "11", "12", "13","14","15"])
        "results/seed2_thresh0.9_h2-0.1_vars0.01/ldsc_enrichment/fdr_heatmap.png"

rule trait_list: 
    input:
       "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/genetic_correlations/geno_correlations.simplified.txt", "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/heritability_estimates/ukbb.tsv",
    output:
        "trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.studies.tsv","trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.names.tsv"
    shell:
        """
        Rscript src/getSelectionList.R --corr_dat {input[0]}  --trait_list {input[1]} --output ./trait_selections/ --nongender_specific --num_samples 1 --start_seed {wildcards.seedn} --threshold {wildcards.thresh} --h2 {wildcards.h2} --names
        """
rule extract_snps: #finds the variants that meet our signal threshold. here just looking at the LDSC ones.
    input:
       trait_list="trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.studies.tsv",
    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.union.txt"
    shell:
        """
            ml python/3.7-anaconda
            python src/unionVariants.py --gwas_list {input.trait_list}  --output {output} --type ldsc --pval {wildcards.pval} --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
        """


rule filter_1KG: #filter those lists for multi-allelic snps, indels, ambiguous snps, etc.
    input:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.union.txt"

    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.ids.txt",
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.1000G.txt"
        
    shell:
        """
        bash src/variant_lookup.sh {input} {output[0]}
        bash src/snp_cleanup.sh {output[0]} {output[1]}
        """

rule prune: #reduce it to a pruned list
    input:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.1000G.txt"
    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/500kb.0.04r2.prune.in"
    
    params: "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/500kb.0.04r2"

    shell:
        """
        plink2 --bfile /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --indep-pairwise 500kb 0.04 --extract {input} --out {params};
        """

rule ids_to_rsids:
    input:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/500kb.0.04r2.prune.in",

    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.pruned_rsids.txt"
    shell: #Only applies if using the LDSC variants, which here we are sadly.    
        """
        bash src/variant_to_rsid.sh {input} {output}
        """

rule extract_sumstats: #get out the z-scores
    input:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.pruned_rsids.txt",        
        "trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.studies.tsv"
    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.se.tsv",
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.n.tsv",
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.beta.tsv"
    params:
        gwas_dir="/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/unzipped", 
        outfile="gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}",
        type="std"
    shell: 
        """
        ml python/3.7-anaconda;
        python src/quickGWASIter.py  --type {params.type} --output {params.outfile} --gwas_list {input[1]} --snp_list {input[0]} --extension ".both_sexes.tsv" --gwas_dir {params.gwas_dir}
        """
rule hapmap_reference: #get the list of hapmap snps for extraction, omitting HLA region
    input:
    output:
        "data/hm3_no_hla.txt"
    shell:
        """
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz -p data/
        tar -xvzf data/weights_hm3_no_hla.tgz
        for i in {{1..22}}; do zcat data/weights_hm3_no_hla/weights.${{i}}.l2.ldscore.gz | tail -n +2 | awk '{{print $1":"$3"\t"$2}}' >> {output}; done
        """

rule hapmap_extract: #Pull the hapmap3 snps from the LDSC summary stats. This step takes a bit more memory, 10 GB at least.
    input:
        "data/hm3_no_hla.txt",
        "trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.studies.tsv"
    output:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/full_hapmap3_snps.z.tsv",
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/full_hapmap3_snps.n.tsv"
    params:
        "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/",
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/full_hapmap3_snps"
    shell: 
        """
        ml python/3.7-anaconda;
        python src/quickGWASIter.py  --type ldsc --output {params[1]} --gwas_list {input[1]} --snp_list {input[0]} --extension ".both_sexes.tsv" --gwas_dir {params[0]}
        """
rule factorize:
    input:
        "gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}.z.tsv", 
        "trait_selections/seed{seedn}_thresh{thresh}_h2-{h2}.{custom}.names.tsv"
    output:
        expand("results/seed{{seedn}}_thresh{{thresh}}_h2-{{h2}}_vars{{pval}}.{{custom}}/factorization/{ot}.{{f_type}}.{{seedn}}{ext}", ot=output_type, ext=[".txt", ".png"])
        #expand("results/seed{{seedn}}_thresh{{thresh}}_h2-{{h2}}_vars{{pval}}/factorization/{ot}.{ft}.{{seedn}}{ext}", ot=output_type, ft = factorization_type, ext=[".txt", ".png"])
    params:
        "results/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}.{custom}/factorization/"
    run:
        shell("Rscript src/MF_wrapper.R")

#project onto all hapmap 3 snps
rule project:
    input:
        factors="factorization_data/{identifier}.factors.txt",
        variants="gwas_extracts/seed{seedn}_thresh{thresh}_h2-{h2}_vars{pval}/full_hapmap3_snps.z.tsv" #fixed, for now.
    output:
        "results/{identifier}/projected_hapmap3_loadings.txt"
    params:
    run:
        shell("Rscript {src_path}/projectSumStats.R --output {output} --factors {input.factors} --sumstats {input.variants} --id_type 'RSID'")

checkpoint prep_enrichment: #format the outputed factors for enrichment analysis
    input:
        hapmap_list="/work-zfs/abattle4/ashton/reference_data/hapmap_chr_ids.txt",
        projections="results/{identifier}/projected_hapmap3_loadings.txt",
        sample_counts="/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/full_hapmap3_snps.n.tsv"

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
        "ldsc_reference/Multi_tissue_gene_expr.ldcts",
        expand("ldsc_reference/weights_hm3_no_hla/weights.{chr}.l2.ldscore.gz", chr = range(1,23))
    shell:
        """
            mkdir -p ldsc_reference
            cd ldsc_reference
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
        ldsc_ref ="ldsc_reference/{tis_ref}.ldcts"
    output:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.cell_type_results.txt",
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue.log"
    params:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{factor}.multi_tissue"
    shell:
        """
        cd ldsc_reference
        python /work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
        --h2-cts ../{input.ss} \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out ../{params} \
        --ref-ld-chr-cts ../{input.ldsc_ref} \
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
        "results/{identifier}/ldsc_enrichment_{tis_ref}/full_heatmap.png", "results/{identifier}/ldsc_enrichment_{tis_ref}/fdr_heatmap.png"
    params:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/"
    shell:
        """
            echo {input}
            Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "fdr_sig" --output {output[1]}
            Rscript {src_path}/visualizeLDSC.R --input_dir {params} --plot_type "horizontal" --output {output[0]}
        """


