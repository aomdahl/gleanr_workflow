# gleanr_workflow: tools for up-and-downstream analysis of gleanr factors
GLEANR is a GWAS matrix factorization tool to estimate sparse latent pleiotropic genetic factors. Factors map traits to a distribution of SNP effects that may capture biological pathways or mechanisms shared by these traits.
This repo contains tools for preprocessing of GWAS data, running GLEANR from the command line, and then evaluating input factors.
It also contains files with specific commands used to performa analysis reported in the main publication:

[**Sparse matrix factorization of GWAS summary statistics robust to sample sharing improves detection and interpretation of factors with diverse genetic architectures**](https://www.biorxiv.org/content/10.1101/2024.11.12.623313v1).


## Installing the GLEANR package
This can be done from the companion repository [here](https://github.com/aomdahl/gleanr/tree/main) 

## Repo structure:
 - `rules`: Snakemake rule files to perform GLEANR analysis, including cleaning and harmonizing input GWAS data, estimating and formatting GLEANR inputs, and analyzing GLEANR factors downstream
 - `src`: Scripts used in pre and post-processing of GLEANR results which need global accessibility across the whole directly. Entirely independent from the GLEANR software package. Also contains helpful scripts for running GLEANR directly from the command line.(`gleaner_run.R`)
 - `manuscript_analyses`: Documentation for analysis run in manuscript, including scripts specific to our analysis of 137 diverse UKBB traits, the comparative analysis of FinnGen vs UKBB, and simulations.
   
## GLEANR method:
This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.
Our approach was inspired by work from Yuan He [here](https://github.com/heyuan7676/ts_eQTLs).

## Running GLEANR
Tutorials and vignettes will be posted within the next month. If you'd like to try `gleanr` in the meantime, use the script `src/gleaner_run.R` to run analysis directly on input matrices of summary statistics.

## Analysis from the paper
- For a detailed outline of the steps and scripts used in the cross-cohort analysis of FinnGen and UKBB GWAS, see the document [here](https://github.com/aomdahl/gleanr_workflow/blob/main/manuscript_analyses/FinnGen_vs_UKBB_analysis/full_analysis.md)
- For a detailed outline of the steps and scripts used in the analysis of 137 diverse phenotypes, see TBD
