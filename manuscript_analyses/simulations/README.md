# Simulations of sample sharing in matrix factorization of GWAS
Our analysis includes two types of simulation to investigate the impact of sample sharing on matrix factorization of GWAS summary statistics.
Our performane and evaluation of these simulations is summarized here:

## Type 1: GWAS-based simulations with specified levels of sample sharing
*Paths are set relative to running from `manuscript_analyses/simulations/1_gwas_based`*

In these simulations, we generated synthetic genotypes based on UKBB MAFs, and synthetic phenotypes modelled after height, weight, and BMI in the UKBBK.

The entire simulation can be run with the script `N30000_multi-overlap_SD-1-1-0.01_h2-0-0_bmi-like.SCALED.Aug2024.sh`:
```bash
cd manuscript_analyses/simulations/1_gwas_based
sbatch N30000_multi-overlap_SD-1-1-0.01_h2-0-0_bmi-like.SCALED.Aug2024.sh 
```

We then evaluate the output factorizations in `../../figure_scripts/fig2_simulations_gleaner_performance.Rmd` to generate figures for the text.

## Type 2: Factor-based simulations with non-random added noise
*Paths are set relative to running from `manuscript_analyses/simulations/2_factor_based`*
Most of the necessary commands are wrapped up in several scripts.

**0) Generate simulated underlying data (factors, sample sharing structures)**

The command for this is wrapped up in the script in *Step 2*, but is called with the script `Rscript src/buildSimDataClean.R`.

**1) Specify the settings to test in simulation**

Specified simulation combinations are given by setting file which specifies the number of samples (N), source of MAF, degree of phenotype overlap, number of repeated samplings to perform, which factroization methods to use, how many K to generate, which V and U to use, which heritabilities to scale to, and any shrinkage on the covariance matrix. See files in `manuscript_analyses/simulations/2_factor_based/setting_files/`

**2) Generate the simulations**
This command:
1) Generates the underlying simulation data (U,V,SE,B,C, etc.)
2) Processes simulation settings files
3) Generates the actual GWAS summary statistics from the simulated data based on setting files
4) Build files containing the commands to do matrix factorization on the simulated data

```bash
cd manuscript_analyses/simulations/2_factor_based
bash run_scripts/launch_all_sims.sh
```

**3) Execute the simulations**
Simulations with no sample sharing:
```bash
sbatch run_scripts/no_overlap_gwas.first_third.sh
sbatch run_scripts/no_overlap_gwas.second_third.sh
sbatch run_scripts/no_overlap_gwas.third_third.sh
```
Simulations with 1-block smaple shairng structure:
```bash
sbatch run_scripts/1b_overlap_gwas.first_third.sh
sbatch run_scripts/1b_overlap_gwas.second_third.sh
sbatch run_scripts/1b_overlap_gwas.third_third.sh
```
Simulations with 2-block smaple shairng structure:
```bash
sbatch run_scripts/2b_overlap_gwas.first_third.sh
sbatch run_scripts/2b_overlap_gwas.second_third.sh
sbatch run_scripts/2b_overlap_gwas.third_third.sh
```
**4) Execute the simulations**

We valuate the output factorizations in `../../figure_scripts/fig2_simulations_gleaner_performance.Rmd` to generate figures for the text.
