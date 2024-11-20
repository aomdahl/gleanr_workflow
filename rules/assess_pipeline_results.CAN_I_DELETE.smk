filelist = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/infertility/fall_2022/sept_2022_list.txt" #get this from the YAML
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



#Run the tissue specific analysis on the data
LDSC_REF_DAT="/data/abattle4/aomdahl1/reference_data/ldsc_reference"

rule all:
	input:
		expand("results/infertility/fall_2022/ldsc_enrichment_Multi_tissue_chromatin/{study_id}.multi_tissue.cell_type_results.txt", study_id=full_study_ids)


rule ldscGWASenrichment: #just run for one, then call on the input.
    input:
        ss="gwas_extracts/{identifier}/{study_id}.sumstats.gz",
        tiss_dir=LDSC_REF_DAT+"/{tis_ref}.ldcts"
    output:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{study_id}.multi_tissue.cell_type_results.txt",
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{study_id}.multi_tissue.log"
    params:
        "results/{identifier}/ldsc_enrichment_{tis_ref}/{study_id}.multi_tissue",
        LDSC_REF_DAT,
        "results/{identifier}/ldsc_enrichment_{tis_ref}/"
    shell:
        """
        mkdir -p {params[2]}
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


