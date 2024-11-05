# GLEANR: GWAS latent embeddings accounting for noise and regularization
GLEANER is a matrix factorization method to estimate sparse latent genetic factors of pleiotropy across GWAS.
This repo contains the `gleanr` R package, in addition to helpful pipeline scripts to implement and use the package.

## Installing GLEANR

## Repo structure:
 - `rules`: Snakemake rule files to perform GLEANR analysis, including cleaning and harmonizing input GWAS data, estimating and formatting GLEANR inputs, and analyzing GLEANR factors downstream
 - `src`: Scripts used in pre and post-processing of GLEANR results. Entirely independent from the GLEANR software package. Also contains helpful scripts for running GLEANR directly from the command line.(`gleaner_run.R`)
 - `gleanr`: contains the R package
## GLEANR method:
This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.
Our approach was inspired by work from Yuan He [here](https://github.com/heyuan7676/ts_eQTLs).
