#!/bin/bash
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
ml anaconda; conda activate renv
#K==41,full thing
#Finngen version first
Rscript src/svd_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv \
--uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv \
--trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt \
--outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/finngen_benchmark_svd_k41 \
-K 41
#PanUKBB next
Rscript src/svd_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2_conservative.beta.tsv \
--uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//ukbb_benchmark_2_conservative.se.tsv \
--trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//pheno.names.txt \
--outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/ukbb_benchmark_svd_k41 \
-K 41

##Do with a different K?
#finngen 19, UKBB 17
#13/12 vs 19/17
#Finngen version first
Rscript src/svd_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv \
--uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv \
--trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt \
--outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/finngen_benchmark_svd_k19 \
-K 19
#PanUKBB next
Rscript src/svd_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2_conservative.beta.tsv \
--uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//ukbb_benchmark_2_conservative.se.tsv \
--trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//pheno.names.txt \
--outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/ukbb_benchmark_svd_k17 \
-K 17

#Generous as possible case- both 17:
Rscript src/svd_run.R --gwas_effects /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv \
--uncertainty /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv \
--trait_names /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt \
--outdir /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/finngen_benchmark_svd_k17 \
-K 17