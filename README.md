# Sparse GWAS matrix factorization
This repo contains an ongoing project to develop a flexible, interpretable, and sparse factorization framework to integrate GWAS data across studies and cohorts. We employ a basic alternating least-squares matrix factoriztion algorithm with sparse priors on learned matrices, while accounting for study uncertainty.

The project is an extension of work by Yuan He [here](https://github.com/heyuan7676/ts_eQTLs), customized for genome-wide association studies and accounting for specific nuances of GWAS data. In particular, this involves modifications to make the method more flexible to summary statistics from various sources across many different phenotypes.

As in the eQTL version of this method, it implements an alternating least squares approach with L1 lasso to learn sparse latent genetic factors.

This work is currently in development, but more updates coming soon!

