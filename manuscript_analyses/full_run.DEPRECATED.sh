#!/bin/bash

#######################################################################
#
#  Script for performing download XT-LDSC of PanUKBB phenotypes
#
#	Ashton Omdahl, jan-feb 2024
#
#######################################################################

#File download info from :/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/trait_selection.Rmd
# Commands to download files generated there

#1) Download files
  cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB
  tail -n +2 download_phenos_flat.sh > tmp && mv tmp download_phenos_flat.sh
  tail -n +2 download_biomarkers_flat.sh > tmp && mv tmp download_biomarkers_flat.sh
  sbatch download_biomarkers_flat.sh
  sbatch download_phenos_flat.sh

#2) Get the EUR relevant information:
	cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/ 
	#zcat full_variant_qc_metrics.txt.bgz |  awk '((NR == 1) || (($11 + 0.0 > 0.9) && ($34 + 0.0 > 0.01) && ($9 == "true"))) {print $1","$2","$3","$4","$5","$6","$14","$20","$26","$32","$38}' | gzip > high_quality_common_variants_EUR.txt.bgz
	#Above version yielded incorrect thing, trying again....
	#zcat full_variant_qc_metrics.txt.bgz |  cut -f 1-6,11,33,34,9 | awk '((NR == 1) || (($7 == "true") && ($8 + 0.0 > 0.9) && ($10 + 0.0 > 0.01))) {print $0}'  | gzip > high_quality_common_variants_EUR.txt.bgz
	#Third time's the charm?
	#awk -f quick_filter.awk <(zcat full_variant_qc_metrics.txt.bgz | cut -f 1-6,11,33,34,9) > intermediate.tmp
	#awk 'NR==FNR{c[$5]++;next} c[$5]<2' intermediate.tmp intermediate.tmp | gzip >  high_quality_common_variants_EUR.txt.bgz	
	#Nope, just did it in R and it worked fine. This is dumb.
#3) Extract via awk the EUR information only from the summary statistics
  bash scripts/build_extract_scripts.sh
  #Launch the commands into 4 scripts		
  bash submit_extractions.sh
  #ALTERNATIVE option- just do it for the 137 traits we want:
  	#ls *.sh all.runs | sed 's/.extract.sh//g' > all.runs
	#LOOKUP=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/panUKBB_complete.studies.tsv
	#while read p; do
	#    if grep -q $p $LOOKUP; then
	#	echo ${p}.extract.sh
	#    fi
	#done < all.runs
	#bash submit_extractions_selected.sh

#4) Manually compile the manifest list. Found in:
  ls /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/panUKBB_complete.studies.tsv
  

#5) Launch the premunging script to get these into a format that the LDSC munger will recognize
	sbatch scripts/premunge_all.sh
	cd /scratch16/abattle4/ashton/snp_network/custom_l1_framework
	#snakemake --snakefile rules/extract_factorize.smk -j 1 gwas_extracts/panUKBB_complete/panUKBB_complete.1000G.txt --configfile config/panUKBB_complete_config.yaml

        #cd gwas_extracts/panUKBB_complete/
	sed -i 's/NEGLOG10_PVAL_EUR/PVAL/g' gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	#It didn't add the signed sum stats,nor did it add SE. look into this. But for now, do it manually
	sed -i 's/--keep-maf/--keep-maf --signed-sumstats beta_EUR,0/g'  gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	sed -i 's/beta_EUR,0/beta_EUR,0 --se se_EUR/g'  gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	#NOTE: the SE call only works with my moded sumstat munger script.

#6) Split the munge commands across individual jobs I can submit separately
ml snakemake
snakemake --snakefile  rules/ldsr_pairwise.smk  -j 1 gwas_extracts/panUKBB_complete/munge_calls/ldsc_line_runs.mean_platelet_volume.sh --configfile config/panUKBB_complete_config.yaml

#7) Submit the jobs to run on slurm job manager:
 snakemake --snakefile rules/ldsr_pairwise.smk gwas_extracts/panUKBB_complete/missingness_report.tsv  -j 5 --configfile config/panUKBB_complete_config.yaml --profile  profiles/rockfish/

#8) Drop traits we are no longer wanting to include:
#update the file for the trait list
	awk -F "-" '(FNR == NR) {arr[$1]; next} !($2 in arr) {print $0}'  ../ldsr_results/panUKBB_complete/omit_traits_activities_environment_nutritional_social.txt panUKBB_complete.studies.tsv  > panUKBB_selected.studies.tsv
	#Update the config file (in vim)Vx
#8)  Generate the pairwise LDSC commands and run them
	#cleanup previous runs:
	rm ldsr_results/panUKBB_complete/rg_ldsr/*.sh
	snakemake --snakefile rules/ldsr_pairwise.smk -j 1 ldsr_results/panUKBB_complete/rg_ldsr/hdl_cholesterol_ldsc.run.sh  --configfile config/panUKBB_complete_config.yaml
	#Modify them to work with the format, LDSC output formatting was irregular and we need to update to EUR.
	for i in  ` ls ldsr_results/panUKBB_complete/rg_ldsr/*.sh`; do

		sed -i 's/-chr//g' $i
		sed -i 's/UKBB.ALL.ldscore\//UKBB.ALL.ldscore\/UKBB.EUR/g' $i
	done

#9) Run cross trait LDSC:
snakemake --snakefile rules/ldsr_pairwise.smk -j 5 ldsr_results/panUKBB_complete//summary_data/gcov_int.tab.csv  --configfile config/panUKBB_complete_config.yaml --profile profiles/rockfish

#10) Manually evaluate the results, select a threshold to narrow down to ~100 studies
	#in R: trait_selection_part2_postLDSC.Rmd
	#Also update the missingness repoort file
	#Modify the missingness report file
	#awk -F "." '(FNR == NR) {arr[$1];next} ($4 in arr) {print $0}' /scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/rg_filtered_0.7_traits.txt  missingness_report.tsv > missingness_report.SUB.tsv
	#mv missingness_report.tsv missingness_report.COMPLETE.tsv && mv missingness_report.SUB.tsv missingness_report.tsv

#11) Create the file for clumping variants based on the frequency of occurence across the set.
	ml anaconda
	conda activate std
	python src/unionVariants.py --gwas_list gwas_extracts/panUKBB_complete/missingness_report.tsv  --output gwas_extracts/panUKBB_complete/panUKBB_complete.union.txt --type ldsc_custom --pval 1e-5 --gwas_dir gwas_extracts/panUKBB_complete/ --output_counts gwas_extracts/panUKBB_complete/panUKBB_complete.union_freq.txt
	#Prep for the clumping by converting into 1000G compatabile ids
	#back filter to ge the values
	cd gwas_extracts/panUKBB_complete/
	awk '(FNR == NR) {arr[$1]=$2;next} ($2 in arr) {print $1"\t"arr[$2]}' panUKBB_complete.union_freq.txt  panUKBB_complete.ids.txt | sed 's/:[ATCG]:[ACTG,]\+\t/\t/g' >  panUKBB_complete.ids_freq.txt

	awk '(FNR == NR) {arr[$1]=$2;next} ($1 in arr) {print $1"\t"arr[$1]}' panUKBB_complete.ids_freq.txt  panUKBB_complete.1000G.txt >  panUKBB_complete.1000G_freqs.txt

	cat <(echo -e "SNP\tP") panUKBB_complete.1000G_freqs.txt  > t && mv t panUKBB_complete.1000G_freqs.txt

#10) Clump using the UKBB reference
	#Run the clumping:
	ml plink/1.90b6.4
	plink --bfile /scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/ref --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --clump  gwas_extracts/panUKBB_complete/panUKBB_complete.1000G_freqs.txt --out gwas_extracts/panUKBB_complete/panUKBB_complete_clumping_250kb_r2_0.1
	#Convert it into a format that snakemake likes:
	awk '{print $3}'   gwas_extracts/panUKBB_complete/panUKBB_complete_clumping_250kb_r2_0.1.clumped > gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.250kb.0.2r2.prune.in

#11) Extract this set of SNPs
	snakemake --snakefile rules/extract_factorize.smk -j 1 gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.beta.tsv --configfile config/panUKBB_complete_config.yaml
#12) Run GLEANER
	sbatch ../../../custom_l1_factorization/run_scripts/GLEANER_panUKBB_6K_snps_2.sh
#13) Project the factors
	sbatch extractProject.sh
		#based on the #SNPs analyzed by FactorGo, I'm not sure I need to do that....
#14) Build the meta-summary stats
	ml snakemake

	snakemake --snakefile  rules/project_assess.smk  -j 1 results/panUKBB_complete_6.1K/loading_ss_files_Meta/ --configfile config/panUKBB_complete_config.yaml

	snakemake --snakefile  rules/project_assess.smk  -j 1 results/panUKBB_complete_6.1K/loading_ss_files_MetaAdj/ --configfile config/panUKBB_complete_config.yaml

#15) Run the full LDSC
	#Adjusted version, which is what we will use
	snakemake --snakefile  rules/project_assess.smk  -j 7 results/panUKBB_complete_6.1K/MetaAdj_ldsc_enrichment_Multi_tissue_chromatin/factor_global_fdr.heatmap.png --configfile config/panUKBB_complete_config.yaml --profile  profiles/rockfish/
	#Standard version:
	snakemake --snakefile  rules/project_assess.smk  -j 5 results/panUKBB_complete_6.1K/Meta_ldsc_enrichment_Multi_tissue_chromatin/factor_global_fdr.heatmap.png --configfile config/panUKBB_complete_config.yaml --profile  profiles/rockfish/ --dry-run

########Alternative version- run it directly on the factors

