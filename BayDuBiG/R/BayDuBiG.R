# ==============================================================================
# BayDuBiG — Bayesian Durbin-based Bi-level Identification of SVGs
#
# A Bayesian bi-level variable-selection model for identifying spatially variable
# genes (SVGs) in spatial transcriptomics data. The group effect follows an OR
# semantics: a gene is flagged as significant if it belongs to at least one
# selected group (pathway).
#
# Main exported functions:
#   - run_BayDuBiG()        Full analysis pipeline (normalization + weights + MCMC + BFDR)
#   - normalize_expression()
#   - normalize_coordinates()
#   - compute_gaussian_weights()
#
# Author: Tianyi Wang
# ==============================================================================

# ----------------------------- Preprocessing -----------------------------------

#' Normalize the gene expression matrix (Anscombe-type transform + centering)
#'
#' @param expression_matrix Numeric matrix with cells in rows and genes in columns.
#' @return The normalized expression matrix (columns centered to mean 0).
#' @export
normalize_expression <- function(expression_matrix) {
  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix (cells x genes)!")
  }
  if (any(is.na(expression_matrix))) {
    warning("NA values found in the expression matrix; replaced with 0")
    expression_matrix[is.na(expression_matrix)] <- 0
  }

  # Estimate the over-dispersion parameter alpha (quadratic mean-variance relation).
  mean_expr <- colMeans(expression_matrix)
  var_expr  <- apply(expression_matrix, 2, stats::var)
  fit <- stats::lm(var_expr ~ I(mean_expr) + I(mean_expr^2))
  alpha_est <- stats::coef(fit)[3]
  if (is.na(alpha_est) || alpha_est <= 1e-12) {
    alpha_est <- 1
  }

  # Anscombe-type variance-stabilizing transform (log form), then center per gene.
  anscombe_transformed <- log(expression_matrix + 1 / alpha_est)
  centered_matrix <- sweep(anscombe_transformed, 2, colMeans(anscombe_transformed), FUN = "-")

  return(centered_matrix)
}

#' Normalize cell coordinates (Z-score standardization)
#'
#' @param cell_coordinates Matrix or data.frame with cells in rows and coordinates
#'   in columns (at least two columns X/Y).
#' @return The normalized coordinate matrix (mean 0, standard deviation 1).
#' @export
normalize_coordinates <- function(cell_coordinates) {
  if (!is.matrix(cell_coordinates) && !is.data.frame(cell_coordinates)) {
    stop("cell_coordinates must be a matrix or data.frame (cells x coordinates)!")
  }
  if (ncol(as.matrix(cell_coordinates)) < 2) {
    stop("cell_coordinates needs at least two columns (X and Y)!")
  }
  coords_mat <- as.matrix(cell_coordinates)
  normalized_coords <- scale(coords_mat, center = TRUE, scale = TRUE)
  rownames(normalized_coords) <- rownames(cell_coordinates)
  colnames(normalized_coords) <- colnames(cell_coordinates)
  return(normalized_coords)
}

#' Compute Gaussian spatial weights (k-nearest neighbors + Gaussian kernel, row-normalized)
#'
#' @param cell_coordinates Normalized coordinate matrix (cells in rows).
#' @param sigma Gaussian kernel bandwidth (default 0.01).
#' @param k Number of nearest neighbors (default 8).
#' @param epsilon Row-sum threshold below which a row is treated as zero (default 1e-2).
#' @return A row-normalized dense Gaussian weight matrix (cells x cells).
#' @export
compute_gaussian_weights <- function(cell_coordinates, sigma = 0.01, k = 8, epsilon = 1e-2) {
  if (!is.matrix(cell_coordinates)) {
    stop("cell_coordinates must be a matrix (call normalize_coordinates first)!")
  }
  if (sigma <= 0) stop("sigma must be positive!")
  if (k < 1 || k >= nrow(cell_coordinates)) {
    stop(paste("k must be between 1 and", nrow(cell_coordinates) - 1, "!"))
  }

  n_cells <- nrow(cell_coordinates)
  W <- matrix(0, nrow = n_cells, ncol = n_cells)
  rownames(W) <- rownames(cell_coordinates)
  colnames(W) <- rownames(cell_coordinates)

  for (i in seq_len(n_cells)) {
    coord_i <- cell_coordinates[i, ]
    distances <- sqrt(rowSums((t(t(cell_coordinates) - coord_i))^2))
    neighbors <- order(distances)[2:(k + 1)]          # exclude self
    W[i, neighbors] <- exp(-distances[neighbors]^2 / (2 * sigma^2))
  }

  row_sums <- rowSums(W)
  non_zero <- row_sums > epsilon
  W[non_zero, ] <- sweep(W[non_zero, ], 1, row_sums[non_zero], FUN = "/")

  return(W)
}

# ----------------------------- Core analysis -----------------------------------

#' Run the full BayDuBiG analysis pipeline
#'
#' @param raw_expression Raw gene expression matrix (cells x genes).
#' @param raw_coordinates Raw cell coordinates (cells x X/Y).
#' @param gene_group_list Named list: each element is an integer vector of the
#'   group (pathway) indices the gene belongs to. Group indices start at 1; use
#'   integer(0) for a gene that belongs to no group.
#' @param X Covariate matrix (cells x covariates), optional. If NULL only the
#'   spatial coordinates are used.
#' @param sigma Gaussian weight bandwidth (default 0.01).
#' @param k Number of nearest neighbors for the Gaussian weights (default 8).
#' @param iter Total number of MCMC iterations (default 2000).
#' @param burn Number of burn-in iterations (default 1000).
#' @param target_bfdr BFDR control target (default 0.05).
#' @param c_candidates Candidate grid of BFDR thresholds c.
#' @param informative_gamma Whether to use the informative gamma prior (default FALSE).
#' @param a_gamma_on,b_gamma_on Gamma prior when the group effect is active.
#' @param a_gamma_off,b_gamma_off Gamma prior when the group effect is inactive.
#' @param verbose Whether to print progress messages (default TRUE).
#'
#' @return A list containing:
#'   \item{tau_gamma_results}{Per-gene PPI (posterior inclusion probability, named vector).}
#'   \item{gene_results}{A data.frame with one row per gene: PPI (tau_gamma), the
#'     three basis log-likelihoods (loglik_original / loglik_periodic /
#'     loglik_exponential), the selected basis, and the SVG flag.}
#'   \item{svg_status}{Logical vector indicating whether each gene is an SVG.}
#'   \item{svg_gene_names}{Names of the genes identified as SVG.}
#'   \item{optimal_c}{The selected BFDR threshold c.}
#'   \item{basis_per_gene}{The covariate basis selected for each gene ("original" / "periodic" / "exponential").}
#'   \item{mcmc}{The raw MCMC return list.}
#' @export
run_BayDuBiG <- function(
    raw_expression,
    raw_coordinates,
    gene_group_list,
    X = NULL,
    sigma = 0.01,
    k = 8,
    iter = 2000,
    burn = 1000,
    target_bfdr = 0.05,
    c_candidates = seq(0.05, 0.9, by = 0.05),
    informative_gamma = FALSE,
    a_gamma_on = 10, b_gamma_on = 1,
    a_gamma_off = 1, b_gamma_off = 10,
    verbose = TRUE) {

  if (verbose) message("============ Preprocessing ============")

  # 1. Normalize expression and coordinates.
  normalized_expr <- normalize_expression(raw_expression)
  normalized_coords <- normalize_coordinates(raw_coordinates)
  X_coord <- normalized_coords[, 1]
  Y_coord <- normalized_coords[, 2]

  # 2. Gaussian spatial weights.
  W_dense <- compute_gaussian_weights(normalized_coords, sigma = sigma, k = k)
  W_sparse <- Matrix::Matrix(W_dense, sparse = TRUE)

  # 3. Build XWX = [X, W·X] (covariates plus their spatial lags) on the fly.
  if (is.null(X)) {
    XWX <- NULL
  } else {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (nrow(X) != nrow(raw_expression)) {
      stop("X must have the same number of rows as the number of cells!")
    }
    WX <- W_dense %*% X
    XWX <- cbind(X, WX)
  }

  if (verbose) message("============ MCMC sampling ============")
  start_time <- Sys.time()

  mcmc_result <- BayDuBiG_mcmc(
    Y = as.matrix(normalized_expr),
    gene_group_list = gene_group_list,
    XWX = XWX,
    X_coord = X_coord,
    Y_coord = Y_coord,
    weights = W_sparse,
    iter = iter,
    burn = burn,
    informative_gamma = informative_gamma,
    a_gamma_on = a_gamma_on, b_gamma_on = b_gamma_on,
    a_gamma_off = a_gamma_off, b_gamma_off = b_gamma_off
  )

  if (verbose) {
    runtime <- difftime(Sys.time(), start_time, units = "secs")
    message(sprintf("MCMC finished in %.1f seconds", as.numeric(runtime)))
  }

  # 4. Per-gene basis selection: for each gene, pick the covariate basis
  #    (original / periodic / exponential) with the largest gene-level mean
  #    log-likelihood, then use that basis's PPI (matches the paper).
  get_pp <- function(k) as.numeric(mcmc_result[[paste0(k, "_result")]]$tau_gamma_upper)
  get_ll <- function(k) as.numeric(mcmc_result[[paste0(k, "_result")]]$gene_log_likelihood_mean)
  ll_original <- get_ll("original")
  ll_periodic <- get_ll("periodic")
  ll_exponential <- get_ll("exponential")
  best_basis <- apply(cbind(ll_original, ll_periodic, ll_exponential), 1, which.max)
  tau_gamma_upper <- ifelse(best_basis == 1, get_pp("original"),
                     ifelse(best_basis == 2, get_pp("periodic"), get_pp("exponential")))
  if (all(is.na(tau_gamma_upper))) {
    stop("tau_gamma_upper (PPI) is empty; cannot determine SVGs.")
  }
  names(tau_gamma_upper) <- colnames(raw_expression)
  basis_per_gene <- c("original", "periodic", "exponential")[best_basis]
  names(basis_per_gene) <- colnames(raw_expression)

  # 5. BFDR control.
  if (verbose) message("============ BFDR control ============")
  calculate_bfdr <- function(ppi, c) {
    one_minus_ppi <- 1 - ppi
    indicator <- as.numeric(one_minus_ppi < c)
    numerator <- sum(one_minus_ppi * indicator)
    denominator <- sum(indicator)
    if (denominator == 0) return(Inf)
    numerator / denominator
  }

  bfdr_values <- vapply(c_candidates, function(c) calculate_bfdr(tau_gamma_upper, c),
                        numeric(1))
  valid_idx <- which(bfdr_values <= target_bfdr)
  if (length(valid_idx) == 0) {
    stop(paste("No value of c satisfies BFDR <=", target_bfdr,
               "; relax target_bfdr or check the data."))
  }
  optimal_c <- c_candidates[min(valid_idx)]

  one_minus_ppi <- 1 - tau_gamma_upper
  svg_status <- one_minus_ppi < optimal_c
  names(svg_status) <- names(tau_gamma_upper)
  svg_gene_names <- names(tau_gamma_upper)[svg_status]

  # Per-gene result table: PPI, the three basis log-likelihoods, the selected
  # basis, and the SVG flag.
  gene_results <- data.frame(
    gene              = colnames(raw_expression),
    tau_gamma         = tau_gamma_upper,
    loglik_original   = ll_original,
    loglik_periodic   = ll_periodic,
    loglik_exponential = ll_exponential,
    selected_basis    = basis_per_gene,
    is_SVG            = svg_status,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message(sprintf("Selected c = %.2f (BFDR <= %.2f)", optimal_c, target_bfdr))
    message(sprintf("Identified %d / %d SVGs", length(svg_gene_names), length(svg_status)))
  }

  list(
    tau_gamma_results = tau_gamma_upper,
    gene_results = gene_results,
    svg_status = svg_status,
    svg_gene_names = svg_gene_names,
    optimal_c = optimal_c,
    target_bfdr = target_bfdr,
    basis_per_gene = basis_per_gene,
    mcmc = mcmc_result
  )
}
