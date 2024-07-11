######
#need to ensure they are aligned, but whatever.
covar.new <- as.matrix(fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/infertility/fall_2022/summary_data/gcov_int.tab.csv"))
covar.new <- as.matrix(Matrix::forceSymmetric(covar.new))
beta.dat <- as.matrix(fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/infertility/fall_2022/infert_fall22.beta.tsv")[,-1])
# From cvCovEst package, with modification

#' Ledoit-Wolf Linear Shrinkage Estimator
#'
#' @description \code{linearShrinkLWEst()} computes an asymptotically optimal
#'  convex combination of the sample covariance matrix and the identity matrix.
#'  This convex combination effectively shrinks the eigenvalues of the sample
#'  covariance matrix towards the identity. This estimator is more accurate
#'  than the sample covariance matrix in high-dimensional settings under fairly
#'  loose assumptions. For more information, consider reviewing the manuscript
#'  by \insertCite{Ledoit2004;textual}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param sample_cov_mat \code{matrix} covariance matrix object to regularize.
#' Previously: "dat A numeric \code{data.frame}, \code{matrix}, or similar object."
#'
#' @importFrom matrixStats sum2
#'
#' @return A \code{matrix} corresponding to the Ledoit-Wolf linear shrinkage
#'  estimate of the covariance matrix.
#'
#' @references
#'  \==Boileau, Philippe, Nima S. Hejazi, Mark J. van der Laan, and Sandrine Dudoit. 2021.
#'  “Cross-Validated Loss-Based Covariance Matrix Estimator Selection in High Dimensions.”
#'   https://arxiv.org/abs/2102.09715.
#'
#' @examples
#' linearShrinkLWEst(dat = mtcars)
#' @export
linearShrinkLWEst <- function(dat, sample_cov_mat) {
  # get the number of variables and observations
  p_n <- ncol(dat)
  n <- nrow(dat)
  
  # compute the sample covariance matrix
  #Change- pass in the covariance matrix directly
  #sample_cov_mat <- coop::covar(dat)
  
  # define the identity
  idn_pn <- diag(p_n)
  
  # estimate the scalers
  dat <- as.matrix(dat)
  m_n <- matrixStats::sum2(sample_cov_mat * idn_pn) / p_n
  d_n_2 <- matrixStats::sum2((sample_cov_mat - m_n * idn_pn)^2) / p_n
  
  #with standard scale:
  scaled.dat <- scale(dat, center = TRUE, scale = FALSE)
  b_bar_n_2 <- apply(scaled.dat, 1,
                     function(x) {
                       matrixStats::sum2((tcrossprod(x) - sample_cov_mat)^2)
                     }
  )
  
  
  b_bar_n_2 <- apply(
    cvCovEst::safeColScale(dat, center = TRUE, scale = FALSE), 1,
    function(x) {
      matrixStats::sum2((tcrossprod(x) - sample_cov_mat)^2)
    }
  )
  
  
  b_bar_n_2 <- 1 / n^2 * 1 / p_n * sum(b_bar_n_2)
  b_n_2 <- min(b_bar_n_2, d_n_2)
  
  # compute the estimator
  estimate <- (b_n_2 / d_n_2) * m_n * idn_pn +
    (d_n_2 - b_n_2) / d_n_2 * sample_cov_mat
  return(estimate)
}

