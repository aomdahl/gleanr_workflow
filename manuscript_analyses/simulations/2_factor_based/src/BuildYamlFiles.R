#!/usr/bin/env Rscript

suppressMessages(library(pacman))
suppressMessages(p_load(magrittr, tidyr, dplyr, ggplot2, data.table,optparse))

BuildYAML <- function(row, seed_init) {
  dir <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/"

  ret.list <- list(
    factors = paste0(dir, paste0(unlist(filter(conv, Name == row$V, Category == "factors") %>% select(Directory, File)), collapse = "/")),
    loadings = paste0(dir, paste0(unlist(filter(conv, Name == row$U, Category == "loadings") %>% select(Directory, File)), collapse = "/")),
    maf = paste0(dir, paste0(unlist(filter(conv, Name == row$MAF, Category == "MAF") %>% select(Directory, File)), collapse = "/")),
    K = row$K,
    iter = row$NITER,
    samp_overlap = gsub("_\\*_", paste0("_", as.character(row$N), "_"),
                        paste0(dir, paste0(unlist(filter(conv, Name == row$N_o, Category == "samp_overlap") %>% select(Directory, File)), collapse = "/"))),
    pheno_corr = paste0(dir, paste0(unlist(filter(conv, Name == row$pheno_overlap, Category == "pheno_corr") %>% select(Directory, File)), collapse = "/")),
    test_methods = row$test_methods,
    noise_scaler = 1,
    bic_param = "sklearn_eBIC",
    init = "V",
    seed_init = seed_init,
    herit_scaler = ifelse(is.na(row$HERIT), "continuous", row$HERIT),
    covar_shrinkage = ifelse(is.na(row$SHRINKAGE), -1, row$SHRINKAGE)
  )

  data.frame(names = names(ret.list), vals = unlist(ret.list))
}

# Argument parsing
option_list <- list(
  make_option(c("-p", "--param_file"), type = "character", help = "Path to the parameter CSV file", metavar = "FILE"),
  make_option(c("-o", "--output_path"), type = "character", help = "Output directory", metavar = "DIR"),
  make_option(c("-c", "--commands"), action = "store_true", default = FALSE,
              help = "Print commands to run simulations", metavar = "FLAG"),
  make_option(c("-s", "--seed"), type = "integer", default = 1, help = "Initial seed value", metavar = "INTEGER")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$param_file) || is.null(opt$output_path)) {
  stop("Both --param_file and --output_path are required.", call. = FALSE)
}

param.setting <- fread(opt$param_file) %>%
  mutate(fnames = paste0(V, "_", U, "_MAF-", MAF, "_N-", N, "_RHO-", pheno_overlap, "_No-", N_o, ".yml"))

if (opt$commands) {
  cat("mkdir -p ", opt$output_path, "\n")
}

seed_init <- opt$seed
for (i in 1:nrow(param.setting)) {
  seed_init <- seed_init + 1
  if (opt$commands) {
    cat(
      "bash src/runSimulation.sh ",
      paste0(opt$output_path, '/', param.setting[i,]$fnames),
      paste0(opt$output_path, gsub(pattern = ".yml", replacement = "", x = param.setting[i,]$fnames), "/"),
      "\n"
    )
  } else {
    message("Currently running ", param.setting[i,]$fnames)
    tab <- BuildYAML(param.setting[i,], as.character(seed_init))
    write.table(
      tab,
      file = paste0(opt$output_path, "/", param.setting[i,]$fnames),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = ","
    )

  }
}
