# load modules
shell.prefix("ml anaconda;") 
# configurations
src_path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
ldsc_ref_dat="/data/abattle4/aomdahl1/reference_data/ldsc_reference"
#Sample counts need to come from somehwere...
#For later:
output_type =["loadings", "factors", "weights"]
#PLEASE put config information into the config file
#config=dict()
#config['file_ext']=".sumstats.gz"
#config['pval_thresh']=1e-4
#make things better with the checkpoint tutorial:https://evodify.com/snakemake-checkpoint-tutorial/
#We assume that you have a list of phenotypes already, and you've converted everything to a friendly LDSC format
#config["pval_thresh"] from yamls files.
pval_thresh=config["pval_thresh"]
ruleorder: hapmap_extract > extract_sumstats

rule all:
	input: "gwas_extracts/saige_benchmark/saige_benchmark.union.1e-5.txt"
       #trait_list="trait_selections/{identifier}/{handle}.studies.tsv",

rule extract_snps: #finds the variants that meet our signal threshold. here just looking at the LDSC ones.
    input:
       trait_list="gwas_extracts/{identifier}/missingness_report.tsv"
       #trait_list=config['filelist'] #this is NOT the right list if you've done the LDSC extraction. its quite different actually
    output:
        "gwas_extracts/{identifier}/{handle}.union.txt"
    shell:
        """
	    conda activate std
            python src/unionVariants.py --gwas_list {input.trait_list}  --output {output} --type ldsc_custom --pval {pval_thresh} --gwas_dir gwas_extracts/{wildcards.identifier}/
            #python {src_path}/unionVariants.py --gwas_list {input.trait_list}  --output {output} --type ldsc --pval config["pval_thresh"] --extension config["file_ext"] --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
        """
#add a quick check in here- how many predominantly zero? just need to warn aboutthis.... 

rule filter_1KG: #filter those lists for multi-allelic snps, indels, ambiguous snps, etc.
    input:
        "gwas_extracts/{identifier}/{handle}.union.txt"

    output:
        "gwas_extracts/{identifier}/{handle}.ids.txt",
        "gwas_extracts/{identifier}/{handle}.1000G.txt"
        
    shell:
        """
	echo "This step may be sub-optimal, but does okay. Recommend reviewing in your specific application"
        bash {src_path}/variant_lookup.sh {input} {output[0]}
        bash {src_path}/snp_cleanup.sh {wildcards.handle} {output[0]} {output[1]}
        """

rule prune: #reduce it to a pruned list
    input:
        "gwas_extracts/{identifier}/{handle}.1000G.txt"
    output:
        "gwas_extracts/{identifier}/{handle}.250kb.0.2r2.prune.in"
    
    params: "gwas_extracts/{identifier}/{handle}.250kb.0.2r2"

    shell:
      """
       ~/.bin/plink-ng/plink2 --bfile /scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/ref --indep-pairwise 250kb 0.2 --extract {input} --out {params};
      """

rule ids_to_rsids:
    input:
        "gwas_extracts/{identifier}/{handle}.250kb.0.2r2.prune.in",
	"gwas_extracts/{identifier}/{handle}.union.txt"

    output:
        "gwas_extracts/{identifier}/{handle}.pruned_rsids.intermediate.txt",
        "gwas_extracts/{identifier}/{handle}.pruned_rsids.txt",
	
    shell: 
        """
        bash {src_path}/variant_to_rsid.sh {input[0]} {output[0]}
        #make sure they are the same length
	LEN_I=`wc -l {input} | cut -f 1 -d " " `
	LEN_O=`wc -l {output} | cut -f 1 -d " " `
	if [ $LEN_O -gt $LEN_I ]
		then
		echo "Warning: different lengths..."
		#exit 1
	fi
	#MAP RSIDs back to the original list:
	conda activate renv
	Rscript src/map_rsid_to_original_list.R {output[0]} {input[1]} {output[1]}
	"""

rule extract_sumstats: #get out the z-scores
    #we also need the summary stats, but not accounting for that here in input
    input:
        "gwas_extracts/{identifier}/{handle}.pruned_rsids.txt",        
        #"trait_selections/{handle}.studies.tsv"
        "gwas_extracts/{identifier}/missingness_report.tsv"
    output:
        "gwas_extracts/{identifier}/{handle}.se.tsv",
        "gwas_extracts/{identifier}/{handle}.n.tsv",
        "gwas_extracts/{identifier}/{handle}.beta.tsv"
    params:
        outfile="gwas_extracts/{identifier}/{handle}",
	gwas_dir="gwas_extracts/{identifier}/",
        type="ldsc_custom"
    shell: 
        """
        conda activate std
	python src/quickGWASIter.py  --type {params.type} --output {params.outfile} --gwas_list {input[1]} --snp_list {input[0]} --extension {config[file_ext]} --gwas_dir {params.gwas_dir}
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

#This is probably overkill- I think everything we need is in the existing LDSC file already, since we subsetted to HM3 at the start.
#can probably do without this.
rule hapmap_extract: #Pull the hapmap3 snps from the LDSC summary stats. This step takes a bit more memory, 10 GB at least.
    input:
        "data/hm3_no_hla.txt",
        #"trait_selections/{identifier}.studies.tsv"
        "gwas_extracts/{identifier}/missingness_report.tsv",
	"gwas_extracts/{identifier}/{handle}.beta.tsv"
    output:
        #"gwas_extracts/{identifier}/full_hapmap3_snps.z.tsv",
        "gwas_extracts/{identifier}/{handle}.full_hapmap3_snps.n.tsv"
    params:
        outfile="gwas_extracts/{identifier}/{handle}.full_hapmap3_snps",
        gwas_dir="gwas_extracts/{identifier}/",
        type="ldsc_custom"
    shell: 
        """
        conda activate std
	    python src/quickGWASIter.py  --type {params.type} --output {params.outfile} --gwas_list {input[1]} --snp_list {input[0]} --extension {config[file_ext]} --gwas_dir {params.gwas_dir} 
        """
        
#this will be the last big change
rule factorize:
    input:
        input_gwas="gwas_extracts/{identifier}/{identifier}.z.tsv", 
        input_names="trait_selections/{identifier}.names.tsv",
        input_se = "gwas_extracts/{identifier}/{identifier}.se.tsv"
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
	#The way to do the scaling correctly is to scale the data....
