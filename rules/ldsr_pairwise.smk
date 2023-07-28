###
#This set of rules is to perform LD-score regression on each GWAS pairwise in order to estimate key statistics (genomic inflation and covariance due to cohort overlap
###


#filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/infertility/fall_2022/sept_2022_list.txt" #get this from the YAML
filelist=config["filelist"]
#filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/saige_benchmark.studies.txt"
#filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/finngen_benchmark.studies.txt"
#filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/udler_original.fall2022.txt"
#filelist="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/ukbb_seed2_fall2022.txt"
shell.prefix("ml anaconda;")
#prep the summary stats for LDSC analysis.

#Split the calls and run separately for each trait so we get the sum stats the way we want it...
def getTraitNames(fl):
    traits = list()
    with open(fl, 'r') as istream:
        for l in istream:
            l = l.strip().split()
            traits.append(l[1])
    return(traits)

def getGWASNames(fl):
    gwas = list()
    import os
    with open(fl, 'r') as istream:
        for l in istream:
            l = l.strip().split()
            basename = os.path.basename(l[0])
            gwas.append(os.path.splitext(basename)[0])
    return(gwas)
trait_list=getTraitNames(filelist)

gwas_list = getGWASNames(filelist)
full_study_ids = [gwas_list[i] + "." + trait_list[i] for i in range(0,len(gwas_list))]


rule all:
	#input:expand("ldsr_results/saige_benchmark/rg_ldsr/{trait_id}_ldsc.run.sh", trait_id = trait_list)
	input: "ldsr_results/finngen_benchmark/summary_data/gcov_int.png"
		#"ldsr_results/saige_benchmark/summary_data/gcov_int.png"
		#expand("ldsr_results/ukbb_GWAS_h2-0.1_rg-0.9/{trait_id}.log", trait_id = trait_list)
		#expand("ldsr_results/saige_benchmark/{trait_id}.log", trait_id = trait_list)
		#expand("ldsr_results/udler_original/{trait_id}.log", trait_id = trait_list)
		#ldsr_results/infertility/fall_2022/
		#expand("ldsr_results/infertility/fall_2022/ldsc_enrichment_Multi_tissue_chromatin/{study_id}.multi_tissue.log", study_id=full_study_ids)

# This creates the bash file to clean up all the summary statistics, based on an input file of GWAS to test on, and output them in a friendly format for LDSC.
#
        #study_list="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/infertility/fall_2022/sept_2022_list.txt", 
        #study_list="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/ukbb_seed2_fall2022.txt", #this file has format filepath pheno_name genome_build sample_count cohort_name

rule mungeCommands:
    input: study_ref = "trait_selections/{identifier}.studies.txt", snp_ref = "/data/abattle4/aomdahl1/reference_data/hm3_snps_ldsc_ukbb.tsv"

    output: "gwas_extracts/{identifier}/munge_sumstats.all.sh"
    params:
        "gwas_extracts/{identifier}/"
    shell:
        """
        conda activate std
        python src/munge_precursor.py --study_list {input.study_ref} --merge_alleles {input.snp_ref} --output {params} --mod --keep_maf
        #we want the mod version...
        """

#This splits those commands up into multiple files run 
rule splitMungeCommands:
    input: "gwas_extracts/{identifier}/munge_sumstats.all.sh"
    output:
        expand("gwas_extracts/{{identifier}}/munge_calls/ldsc_line_runs.{runid}.sh", runid=trait_list)
    params:
        "gwas_extracts/{identifier}/munge_calls/ldsc_line_runs."
    run:
        import re
        ex_list = list()
        print(str(input[0]))
        with open(str(input[0]), 'r') as istream:
            for l in istream:
                txt = l.strip()
		try:
                    tag = re.search("--out [-\w\/\.-]+\.([\w_-]+) --", txt)
                    with open(params[0] + str(tag.group(1)) + ".sh", 'w') as ostream:
                        ostream.write(txt)
			#print("wrote out", txt)
			print("to", (params[0] + str(tag.group(1)) + ".sh"))	
                        if "INTERMEDIATE." in txt:
                            intermediate_file=re.search("--sumstats (.+) --out", txt)
                            ostream.write("\n#rm " + intermediate_file.group(1))
                    ex_list.append(tag)
                except:  
                    print(txt)
                    print("Error with REGEX- debug")


#use the expand to get through the entire list?
#th4 names basically come from the full file name less the extension. Not sure this is the way we want to do it...
#NEed to test this out...
#This actually converts the summary stats into an LDSC-friendly format..
rule ldscMunge:
    input: "gwas_extracts/{identifier}/munge_calls/ldsc_line_runs.{trait_id}.sh"
    output: 
        "gwas_extracts/{identifier}/{gwas_id}.{trait_id}.sumstats.gz"
    shell:
	    """
        conda activate py27
	    bash {input}
	    """

gwas_list = getGWASNames(filelist)
full_study_ids = [gwas_list[i] + "." + trait_list[i] for i in range(0,len(gwas_list))]

#Do a check on the quality of the summary statistics- how much missingness do they each have?
rule checkMissingness:
	input: 
		studies=expand("gwas_extracts/{{identifier}}/{study_id}.sumstats.gz",study_id=full_study_ids),
		order="trait_selections/{identifier}.studies.txt"
	output: "gwas_extracts/{identifier}/missingness_report.tsv", "gwas_extracts/{identifier}/sample_sd_report.tsv"
	params:
		"gwas_extracts/{identifier}/"
	shell:
		"""
		conda activate renv; 
		Rscript src/qc_check_munged_ss.R --gwas_dir {params[0]} --output {params[0]} --gwas_ext {config[file_ext]} --missing_thresh 0.1 --trait_order {input.order}
		"""
		


#This lines up all the files we wish to compute pairwise covariance for
rule generateLDSCCommands:
    input: expand("gwas_extracts/{{identifier}}/{study_id}.sumstats.gz",study_id=full_study_ids)
    output:
        traitfile = "ldsr_results/{identifier}/pairwise.traits.txt",
    
    run:
        #Make a file with all the names for running
        import os
        wd=os.getcwd()
        with open(output.traitfile, 'w') as ostream:
            for i in input:
                ostream.write(wd + "/" + i + '\n') 
        #Generate all the commands. Puts each one in its own file...
#This creates and calls the commands
rule pairwiseCommands:
    input:
        "ldsr_results/{identifier}/pairwise.traits.txt","gwas_extracts/{identifier}/missingness_report.tsv"
    output:
        expand("ldsr_results/{{identifier}}/rg_ldsr/{trait_id}_ldsc.run.sh", trait_id = trait_list)
    shell:
	    """
	     fp=`pwd` 
	     mkdir -p $fp/ldsr_results/{wildcards.identifier}/rg_ldsr/
	     bash src/all_pairwise_r2g.sh {input[0]}  $fp/ldsr_results/{wildcards.identifier}/rg_ldsr/
	    """

#This actually runs the LDSC that is needed. Allows us to just run single ones as needed.
rule pairwiseLDSC:
    input:
        "ldsr_results/{identifier}/rg_ldsr/{trait_id}_ldsc.run.sh"
    output:
        "ldsr_results/{identifier}/rg_ldsr/{trait_id}.log"
    shell:
        """
		d=`pwd`
       		bash {input}
		cd $d
        """

rule tabularizePairwiseLDSC:
    input:
        "ldsr_results/{identifier}/rg_ldsr/{trait_id}.log"
    output:
        "ldsr_results/{identifier}/rg_ldsr/tabular/{trait_id}.ldsc_report.csv",
        "ldsr_results/{identifier}/rg_ldsr/tabular/{trait_id}.rg_report.tsv"

    params:
        out="ldsr_results/{identifier}/rg_ldsr/tabular/"
    shell:
        """
            conda activate std
            echo {input}
            python src/parseLDSCr2g.py --input_file {input} --output {params.out}
            grep -A 1000 "gcov_int_se" {input} | grep -v "Analysis finished" | grep -v "Total time" > {params.out}{wildcards.trait_id}.rg_report.tsv
        """


rule synthesizesPairwiseLDSC:
    input:
        expand("ldsr_results/{{identifier}}/rg_ldsr/tabular/{trait_id}.ldsc_report.csv", trait_id =trait_list),
        expand("ldsr_results/{{identifier}}/rg_ldsr/tabular/{trait_id}.rg_report.tsv", trait_id = trait_list)
    output:
        "ldsr_results/{identifier}/summary_data/gcov_int.png",
        "ldsr_results/{identifier}/summary_data/gcov_int.tab.csv"
    params:
        input="ldsr_results/{identifier}/rg_ldsr/tabular/",
        out="ldsr_results/{identifier}/summary_data/"
    shell:
        """
        conda activate renv
	Rscript src/processLDSCRgResults.R --input_path {params.input} --which_data "ALL" --output {params.out} --filter_se
	conda deactivate
	#TODO: make sure all the fdr corrections and Se filtering settings line up for all tests, make sense
        """

#This rule is a helper for jumping into the next stream of snakemake commands.
#Still unsure if this is the best way to go though....
#Summarize the distribution of the LDSC outputs- is there signal there? Is there missingness?
#This will determine if you need to extract from the files directly, or can use these hm3 simplified ones.
