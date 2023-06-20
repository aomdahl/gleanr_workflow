alphas=["0.1", "0.2","0.3","0.4", "0.5", "0.6","0.7"]
lambdas=["0.001", "0.01", "0.05", "0.1", "0.2"]
shell.prefix("ml anaconda; conda activate renv;")
rule init_k:
    input:
        gwas="gwas_extracts/{identifier}/B.tsv",
        se = "gwas_extracts/{identifier}/SE.tsv",
	lambda_gc="gwas_extracts/{identifier}/Lambda_gc.tsv"
	
    output:
        "results/{identifier}/factorization/init_k.txt"

    params:
        odir="results/{identifier}/factorization/"
    shell:
        """
        Rscript src/MF_wrapper.R --gwas_effects {input.gwas} --uncertainty {input.se} --weighting_scheme "B_SE" \
        --alphas "0.1,0.2,0.3,0.5,0.6,0.7" --lambdas "0.001,0.01,0.05,0.1,0.2" --scaled_sparsity --output {params.odir} --fixed_first \
         --calibrate_k --genomic_correction {input.lambda_gc}
        """

search_iter=8
rule factorize_explore_broad:
    input:
    	gwas="gwas_extracts/{identifier}/B.tsv",
        se = "gwas_extracts/{identifier}/SE.tsv",
        lambda_gc="gwas_extracts/{identifier}/Lambda_gc.tsv"
    output:
        expand("results/{{identifier}}/factorization/K{{k}}_A{alpha}_L{lambdas}_B_SE_{iter}.png", alpha=alphas, lambdas = lambdas, iter = list(range(1,search_iter+1))),
        "results/{identifier}/factorization/K{k}_performance_summary_file.txt"
    params:
        odir="results/{identifier}/factorization/K{k}_"
    shell:
        """
	#this is 30 x 10 exploration points. So something on the order of 15-20 hoursa run?
	#add the cores function....
        Rscript src/MF_wrapper.R --gwas_effects {input.gwas} --uncertainty {input.se} --weighting_scheme "B_SE" \
        --alphas "0.05,0.15,0.3,0.45,0.6" --lambdas "0.001,0.01,0.05,0.1,0.3" --scaled_sparsity --output {params.odir} --fixed_first \
        --overview_plots --subsample 2500 --niter 15 --nfactors {wildcards.k} --genomic_correction {input.lambda_gc} --parameter_optimization {search_iter} --cores 1
        """

rule set_k_looks:
    input: "results/{identifier}/factorization/init_k.txt"
    output: "results/{identifier}.out"
    run:
        with open(input, 'r') as istream:
            k = int(istream.readline().strip())
        if k > 10:
            k_opts = [k - 10, k-5, k, k+5, k+10]
        elif k > 5:
            k_opts = [k-5, k, k+5, k+10]
        elif k < 5:
            k_opts = [k-2, k, k+2, k+4]
        else:
            k_opts = [5,10,15,20,25,30]

k_opts = [3,8,13,18]
#This rule will generate the commands to run a maximal run.
rule factorize_max:
    input:
        expand("results/{{identifier}}/factorization/K{k}_performance_summary_file.txt", k=k_opts)
    output:
        "results/{identifier}/factorization/max_summary_file_iter1.txt"

    run:
        import pandas as pd
        tab_list = []
        for kinfo in input:
            with open(kinfo, 'r') as istream:
                tab_list.append(pd.read_table(kinfo))

        data = pd.concat(tab_list)
        data.to_csv(output)
       	#script("Rscript src/...")

            
#This rule is a backup option if doing manually and something goes wrong....
rule max_performance:
	input:
		expand("results/{{identifier}}/factorization/K{k}_A{alpha}_L{lambdas}_B_SE_{iter}.png", k=k_opts, alpha=alphas, lambdas = lambdas, iter = list(range(1,search_iter+1)))	
	output:
		"results/{identifier}/factorization/overall_performance_summary.txt"
	params:
		path = "results/{identifier}/factorization/"	
	shell:
		"""
			ls {params.path}/*.factors.txt | sed 's/_B_SE.*.factors.txt//g' | sort | uniq  > {params.path}/global_performance_search.txt 
			Rscript src/evaluateRunsMany.R --input_path {params.path} --run_list {params.path}/global_performance_search.txt --output {output}
		"""

