# gleanr_workflow: tools for up-and-downstream analysis of gleanr factors
This repo contains tools for preprocessing of GWAS data, code for performing GWAS matrix factorization via GLEANR from the command line, and then evaluating input factors using various enrichment tests and visualizations.

[GLEANR](https://github.com/aomdahl/gleanr/tree/main) is a GWAS matrix factorization tool to estimate sparse latent pleiotropic genetic factors. Factors map traits to a distribution of SNP effects that may capture biological pathways or mechanisms shared by these traits. This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty, taking inspiration from generalized least squares by transforming regression terms to account for spurious covariance structure while simultaneously applying lasso regularization. (For a more complete description of the method, see our [publication]( https://www.cell.com/ajhg/fulltext/S0002-9297\(25\)00273-3).)

It also contains files with specific commands used to perform the analysis reported in the main publication:

[**Sparse matrix factorization robust to sample sharing across GWASs reveals interpretable genetic components**](https://www.cell.com/ajhg/fulltext/S0002-9297\(25\)00273-3)

## Installing the GLEANR package
This can be done from the companion repository [here](https://github.com/aomdahl/gleanr/tree/main) 

## Repo structure:
 - `rules`: Snakemake rule files to perform GLEANR analysis, including cleaning and harmonizing input GWAS data, estimating and formatting GLEANR inputs, and analyzing GLEANR factors downstream
 - `src`: Scripts used in pre and post-processing of GLEANR results which need global accessibility across the whole directly. Entirely independent from the GLEANR software package. Also contains helpful scripts for running GLEANR directly from the command line.([gleanr_run.R](https://github.com/aomdahl/gleanr_workflow/blob/main/src/gleanr_run.R))
 - `manuscript_analyses`: Documentation for analysis run in manuscript, including scripts specific to our analysis of 137 diverse UKBB traits, the comparative analysis of FinnGen vs UKBB, and simulations.

*Note that currently, the `Snakemake` file in the top directory is not used but is simply a placeholder. The key working functionality of this pipeline resides in the Snakemake rules files in `rules`*
## GWAS factor analysis workflow overview:
In order to both run and interpret GWAS factor analysis,we provide a general workflow and code for performing the following tasks:
1) Harmonizing, munging, and cleaning up input GWAS summary statistics of various formats into a standardized, LDSC-like format for downstream processing
   * This takes as input a list of GWAS summary statistics files (and their paths), and is described by the `mungeCommands`, `splitMungeCommands`, `ldscMunge` and `checkMissingness` rules in [ldsr_pairwise.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/ldsr_pairwise.smk)
2) Estimating covariance inflation due to sample sharing via pairwise XT-LDSC of all input GWAS:
   * This is described by the `generateLDSCCommands`, `pairwiseCommands`, `pairwiseLDSC`, `tabularizePairwiseLDSC`, and `synthesizesPairwiseLDSC` rules in [ldsr_pairwise.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/ldsr_pairwise.smk)
   * The output are tables and matrices containing summaries of the LDSC statistics from pairwise analysis of input GWAS, some of which are direct inputs to `gleanr`.
3) Extracting the relevant GWAS data for facorization analysis.
   * This is described in the rules of [extract_factorize.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/extract_factorize.smk)
4) Performing matrix factorization via `gleanr`
   This can be done using the `gleanr_run.R` script (see the **Running GLEANR** section below)
5) Interpretation of estimated factors, such as by
   * Test factors for tissue-specific marker enrichment (described in [tissue_enrichment.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/tissue_enrichment.smk))
   * Prioritizing top SNPs and their genes based on the distribution of factor effects  (described in [tissue_enrichment.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/tissue_enrichment.smk), rules `top_snps_by_factor` and `top_genes_by_factor`
   * Test factors for gene set enrichment (described in [tissue_enrichment.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/tissue_enrichment.smk), rules `gene_set_enrichment`, `join_enrichment_dat`, `complete_enrichment_analysis`.
   * Characterize factor genetic architecture, described by rule `selective_pressure` in  [interpret_factors.smk](https://github.com/aomdahl/gleanr_workflow/blob/main/rules/interpret_factors.smk)

## Running GLEANR
A vignette for using the R package directly has been posted on the [gleanr repository](https://github.com/aomdahl/gleanr/blob/main/vignettes/gleanr-basic.Rmd).

If you'd like to run `gleanr` directly from the command line, use the script `src/gleaner_run.R` to run analysis directly on input matrices of summary statistics. This can be done using the same data as used in the posted vignette, and available [here](https://github.com/aomdahl/gleanr/tree/main/inst/extdata):
```
SIM_TEST=/path/to/sim/data
Rscript src/gleaner_run.R --gwas_effects ${SIM_TEST}/sim1.effect_sizes.txt --uncertainty ${SIM_TEST}/sim1.std_error.txt --trait_names ${SIM_TEST}/sim1.pheno_names.txt --outdir tutorial_test/test1 \
 --fixed_first -K "GRID" -v 1 \
 --covar_matrix ${SIM_TEST}/sim1.c_matrix.txt \
 --WLgamma Strimmer \
 --covar_se_matrix ${SIM_TEST}/sim1.c_se_matrix.txt --converged_obj_change 0.005
```
This should yield outputs nearly identical to those from the vignette (with the major difference of here, SNPs are ordered lexiciographically as opposed to ordered as given). For a full explantion of arguments, run `Rscript src/gleanr_run.R --help`. A few highlights:

  * `--gwas_effects`: a tabular file containing GWAS effect size estimates ($\beta$s, $N$ SNPs by $M$ traits, plus a column of SNP IDs and a row of trait names)
  * `--uncertainty`: a corresponding input matrix of GWAS standard error estimates ($\beta$s, $N$ SNPs by $M$ traits, plus a column of SNP IDs and a row of trait names)
  * `--trait_names`: a file with trait names (corresponding to the order of $M$) on separate lines
  * `--outdir`: target destination for output files
  * `--fixed_first`: boolean argument, specifies that first column of V will be unregularized
  * `-K`: Number of factors $K$ to initialize at. `GRID` indicates a grid search will be conducted to minimize BIC 
  * `--covar_matrix`: a tabular file with estimates of correlation due to sample sharing for $M$ input traits. Rows and columns are in the same order as `trait_names`, $M \times M$.
  * `--WLgamma`: specifies the type of shrinkage to apply to the covariance matrix. Default is "STRIMMER", which uses the `covar_se_matrix` to apply shrinkage based on estimate uncertainty. May also be a proportion (between 0 and 1, where 1 implies complete shrinkage).
  * `--covar_se_matrix`: a tabular file corresponding to `covar_matrix` with standard errors of the estimates of correlation due to sample sharing for $M$ input traits.
  * `--converged_obj_change`: Specify the proportion of change in the objective to declare convergence. A value of 0.01 means that convergence occurse when the objective changes by 1% or less from the previous iteration.
    
## Analysis from the paper
For a detailed outline of particular analyses and results from the [publication](https://www.cell.com/ajhg/fulltext/S0002-9297(25)00273-3) (including the relevant scripts), please see the following:
 - Simulation of GWAS summary statistics and benchmarking against existing factorization methods [here](https://github.com/aomdahl/gleanr_workflow/tree/main/manuscript_analyses/simulations)
 - Cross-cohort analysis of FinnGen and UKBB GWAS, see the document [here](https://github.com/aomdahl/gleanr_workflow/blob/main/manuscript_analyses/FinnGen_vs_UKBB_analysis/)
 - Analysis of 137 diverse phenotypes, see [here](https://github.com/aomdahl/gleanr_workflow/tree/main/manuscript_analyses/PanUKBB_analysis)
   
For access to several of the outputs from the published work, including GLEANR factorization results and enrichment findings, see our Zenodo link [here](https://zenodo.org/records/15190916)
