# load modules
shell.prefix("module load gcc/5.5.0; module load R;") 
# configurations
configfile: "config/config.yaml"
src_path="/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/"
#For later:
output_type =["loadings", "factors", "weights"]
#make things better with the checkpoint tutorial:https://evodify.com/snakemake-checkpoint-tutorial/
rule all:
    input:
       # expand("results/seed2_thresh0.9_h2-0.1_vars0.01/ldsc_enrichment/F{factor}.multi_tissue.cell_type_results.txt", factor=["1","2","3","4","5","6","7", "8", "9", "10", "11", "12", "13","14","15"])
        "results/seed2_thresh0.9_h2-0.1_vars0.01/ldsc_enrichment/fdr_heatmap.png"

rule trait_list: 
    input:
       "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/genetic_correlations/geno_correlations.simplified.txt", "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/heritability_estimates/ukbb.tsv",
    output:
        "trait_selections/{identifier}.studies.tsv","trait_selections/{identifier}.names.tsv"
    shell:
        """
        Rscript {src_path}/getSelectionList.R --corr_dat {input[0]}  --trait_list {input[1]} --output ./trait_selections/ --nongender_specific --num_samples 1 --start_seed 2 --threshold 0.9 --h2 0.1 --names
        """
rule extract_snps: #finds the variants that meet our signal threshold. here just looking at the LDSC ones.
    input:
       trait_list="trait_selections/{identifier}.studies.tsv",
    output:
        "gwas_extracts/{identifier}/{identifier}.union.txt"
    shell:
        """
            ml python/3.7-anaconda
            python {src_path}/unionVariants.py --gwas_list {input.trait_list}  --output {output} --type ldsc --pval 1e-5 --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
        """


rule filter_1KG: #filter those lists for multi-allelic snps, indels, ambiguous snps, etc.
    input:
        "gwas_extracts/{identifier}/{identifier}.union.txt"

    output:
        "gwas_extracts/{identifier}/{identifier}.ids.txt",
        "gwas_extracts/{identifier}/{identifier}.1000G.txt"
        
    shell:
        """
        bash {src_path}/variant_lookup.sh {input} {output[0]}
        bash {src_path}/snp_cleanup.sh {output[0]} {output[1]}
        """

rule prune: #reduce it to a pruned list
    input:
        "gwas_extracts/{identifier}/{identifier}.1000G.txt"
    output:
        "gwas_extracts/{identifier}/500kb.0.04r2.prune.in"
    
    params: "gwas_extracts/{identifier}/500kb.0.04r2"

    shell:
        """
        plink2 --bfile /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --indep-pairwise 500kb 0.04 --extract {input} --out {params};
        """

rule ids_to_rsids:
    input:
        "gwas_extracts/{identifier}/500kb.0.04r2.prune.in",

    output:
        "gwas_extracts/{identifier}/{identifier}.pruned_rsids.txt"
    shell: #Only applies if using the LDSC variants, which here we are sadly.    
        """
        bash {src_path}/variant_to_rsid.sh {input} {output}
        """

rule extract_sumstats: #get out the z-scores
    input:
        "gwas_extracts/{identifier}/{identifier}.pruned_rsids.txt",        
        "trait_selections/{identifier}.studies.tsv"
    output:
        "gwas_extracts/{identifier}/{identifier}.se.tsv",
        "gwas_extracts/{identifier}/{identifier}.n.tsv",
        "gwas_extracts/{identifier}/{identifier}.beta.tsv"
    params:
        gwas_dir="/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/unzipped", 
        outfile="gwas_extracts/{identifier}/{identifier}",
        type="std"
    shell: 
        """
        ml python/3.7-anaconda;
        python {src_path}/quickGWASIter.py  --type {params.type} --output {params.outfile} --gwas_list {input[1]} --snp_list {input[0]} --extension ".both_sexes.tsv" --gwas_dir {params.gwas_dir}
        """
rule hapmap_reference: #get the list of hapmap snps for extraction, omitting HLA region
    input:
    output:
        "data/hm3_no_hla.txt"
    shell:
        """
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz -p data/
        tar -xvzf data/weights_hm3_no_hla.tgz
        for i in "factorization_data/{identifier}.factors.txt"{{1..22}}; do zcat data/weights_hm3_no_hla/weights.${{i}}.l2.ldscore.gz | tail -n +2 | awk '{{print $1":"$3"\t"$2}}' >> {output}; done
        """

rule hapmap_extract: #Pull the hapmap3 snps from the LDSC summary stats. This step takes a bit more memory, 10 GB at least.
    input:
        "data/hm3_no_hla.txt",
        "trait_selections/{identifier}.studies.tsv"
    output:
        "gwas_extracts/{identifier}/full_hapmap3_snps.z.tsv",
        "gwas_extracts/{identifier}/full_hapmap3_snps.n.tsv"
    params:
        "/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/",
        "gwas_extracts/{identifier}/full_hapmap3_snps"
    shell: 
        """
        ml python/3.7-anaconda;
        python {src_path}/quickGWASIter.py  --type ldsc --output {params[1]} --gwas_list {input[1]} --snp_list {input[0]} --extension ".both_sexes.tsv" --gwas_dir {params[0]}
        """
rule factorize:
    input:
        input.gwas="gwas_extracts/{identifier}/{identifier}.z.tsv", 
        input.names="trait_selections/{identifier}.names.tsv",
        input.se = "gwas_extracts/{identifier}/{identifier}.se.tsv"
    output:
        "results/{identifier}/factorization/objective.png", "results/{identifier}/factorization/runDat.RData"
    params:
        odir="results/{identifier}/factorization/"
    run:
        """
        Rscript src/MF_wrapper --gwas_effects {input.gwas} --uncertainty {input.se} --weighting_scheme "B_SE" \
        --alphas "0.001,0.01,0.05,0.1,0.2" --lambdas "0.001,0.01,0.05,0.1,0.2" --scaled_sparsity --output {params.odir} --fixed_first \
        --overview_plots --nfactors 15 --niter 50 --trait_names {input.names} 
        """
