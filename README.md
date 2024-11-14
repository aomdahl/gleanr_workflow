# GLEANR: GWAS latent embeddings accounting for noise and regularization
<img align="left" src="gleanr/gleanr_logo.png" width="125">GLEANER is a GWAS matrix factorization tool to estimate sparse latent pleiotropic genetic factors. Factors map traits to a distribution of SNP effects that may capture biological pathways or mechanisms shared by these traits.
This repo contains the `gleanr` R package (in development), in addition to helpful pipeline scripts to implement and use the package.
The bioRxiv preprint describing the `gleanr` method in detail is avaialable [here](https://www.biorxiv.org/content/10.1101/2024.11.12.623313v1).


## Installing GLEANR
This can be done directly from github using the  `devtools` package as follows:
```
devtools::install_github("aomdahl/gleanr/gleanr")
```

## Repo structure:
 - `rules`: Snakemake rule files to perform GLEANR analysis, including cleaning and harmonizing input GWAS data, estimating and formatting GLEANR inputs, and analyzing GLEANR factors downstream
 - `src`: Scripts used in pre and post-processing of GLEANR results. Entirely independent from the GLEANR software package. Also contains helpful scripts for running GLEANR directly from the command line.(`gleaner_run.R`)
 - `gleanr`: contains the R package
## GLEANR method:
This is an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.
Our approach was inspired by work from Yuan He [here](https://github.com/heyuan7676/ts_eQTLs).

## Running GLEANR
Tutorials and vignettes will be posted within the next month. If you'd like to try `gleanr` in the meantime, use the script `src/gleaner_run.R` to run analysis directly on input matrices of summary statistics.
