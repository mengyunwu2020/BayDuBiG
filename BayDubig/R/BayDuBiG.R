# ==============================================================================
# Baydubig Analysis Framework (Complete Pipeline with Normalization)
# Includes: Expression normalization, coordinate scaling, Gaussian weight calculation, MCMC analysis + BFDR control
# Author: Tianyi Wang
# Date: 2025-12-24
# ==============================================================================

# Load required packages (ensure dependencies are complete)
required_packages <- c("Rcpp", "RcppDist", "RcppEigen", "Matrix", "MASS", "DirichletReg", "gtools", "here", "qvalue", "ggplot2")
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

# Compile C++ source code (configure here root directory in advance)
compile_baydubig_cpp <- function(cpp_path = here("Baydubig", "src", "Baydubig_mcmc.cpp")) {
  if (!file.exists(cpp_path)) {
    stop(paste("C++ file not found:", cpp_path))
  }
  Rcpp::sourceCpp(cpp_path, verbose = TRUE)
  message("✅ Baydubig_mcmc compiled successfully")
}

# -------------------------- Preprocessing Functions --------------------------
#' Normalize Gene Expression Matrix (Anscombe Transform + Centering)
#' 
#' @param expression_matrix Numeric matrix of gene expression (rows = cells, columns = genes)
#' @return Normalized expression matrix (log-transformed + centered)
normalize_expression <- function(expression_matrix) {
  # Validate input
  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix (cells × genes)!")
  }
  if (any(is.na(expression_matrix))) {
    warning("NA values detected in expression matrix - replaced with 0")
    expression_matrix[is.na(expression_matrix)] <- 0
  }
  
  message("📈 Starting expression normalization (Anscombe transform + centering)...")
  
  # Step 1: Compute total transcript counts per cell
  total_counts <- rowSums(expression_matrix)
  
  # Step 2: Estimate overdispersion (mean-variance relationship)
  mean_expr <- colMeans(expression_matrix)  # Per-gene mean (columns = genes)
  var_expr <- apply(expression_matrix, 2, var)  # Per-gene variance
  fit <- lm(var_expr ~ I(mean_expr) + I(mean_expr^2))
  alpha_est <- coef(fit)[3]
  
  # Handle near-zero overdispersion
  if (alpha_est <= 1e-12) {
    warning("Overdispersion estimate ≤ 1e-12 - reset to 1")
    alpha_est <- 1
  }
  
  # Step 3: Apply Anscombe's transform for variance stabilization
  anscombe_transformed <- log(expression_matrix + 1 / alpha_est)
  
  # Step 4: Center matrix (per-gene mean = 0)
  centered_matrix <- sweep(anscombe_transformed, 2, colMeans(anscombe_transformed), FUN = "-")
  
  message("✅ Expression normalization completed (overdispersion = ", round(alpha_est, 4), ")")
  return(centered_matrix)
}

#' Normalize Cell Coordinates (Z-score Standardization)
#' 
#' @param cell_coordinates Numeric matrix/data frame of cell coordinates (rows = cells, columns = X/Y)
#' @return Z-score normalized coordinate matrix (mean = 0, SD = 1)
normalize_coordinates <- function(cell_coordinates) {
  # Validate input
  if (!is.matrix(cell_coordinates) && !is.data.frame(cell_coordinates)) {
    stop("cell_coordinates must be a matrix or data frame (cells × coordinates)!")
  }
  if (ncol(cell_coordinates) < 2) {
    stop("cell_coordinates must have at least 2 columns (X and Y)!")
  }
  if (any(is.na(cell_coordinates))) {
    warning("NA values detected in coordinates - replaced with column mean")
    cell_coordinates[is.na(cell_coordinates)] <- colMeans(cell_coordinates, na.rm = TRUE)
  }
  
  message("🌍 Starting coordinate normalization (Z-score standardization)...")
  
  # Convert to matrix and apply Z-score
  cell_coords_mat <- as.matrix(cell_coordinates)
  normalized_coords <- scale(cell_coords_mat, center = TRUE, scale = TRUE)
  
  # Restore row names
  rownames(normalized_coords) <- rownames(cell_coordinates)
  colnames(normalized_coords) <- colnames(cell_coordinates)
  
  message("✅ Coordinate normalization completed (mean = 0, SD = 1)")
  return(normalized_coords)
}

#' Compute Gaussian Spatial Weights (k-nearest neighbors)
#' 
#' @param cell_coordinates Normalized cell coordinate matrix (rows = cells, columns = X/Y)
#' @param sigma Bandwidth parameter for Gaussian kernel (default = 0.01)
#' @param k Number of nearest neighbors (default = 5)
#' @param epsilon Small value to avoid division by zero (default = 1e-2)
#' @return Dense Gaussian weight matrix (rows/columns = cells, row-normalized to sum to 1)
compute_gaussian_weights <- function(cell_coordinates, sigma = 0.01, k = 5, epsilon = 1e-2) {
  # Validate input
  if (!is.matrix(cell_coordinates)) {
    stop("cell_coordinates must be a matrix (use normalize_coordinates first)!")
  }
  if (sigma <= 0) {
    stop("sigma must be a positive number (Gaussian bandwidth)!")
  }
  if (k < 1 || k >= nrow(cell_coordinates)) {
    stop(paste("k must be between 1 and", nrow(cell_coordinates) - 1, "(number of cells)!"))
  }
  
  message("⚖️ Computing Gaussian spatial weights (k = ", k, ", sigma = ", sigma, ")...")
  
  # Number of cells
  n_cells <- nrow(cell_coordinates)
  
  # Initialize dense weight matrix
  W_gaussian_dense <- matrix(0, nrow = n_cells, ncol = n_cells)
  rownames(W_gaussian_dense) <- rownames(cell_coordinates)
  colnames(W_gaussian_dense) <- rownames(cell_coordinates)
  
  # Calculate weights for each cell
  for (i in 1:n_cells) {
    # Coordinates of current cell
    coord_i <- cell_coordinates[i, ]
    
    # Euclidean distance to all other cells
    distances <- sqrt(rowSums((t(t(cell_coordinates) - coord_i))^2))
    
    # Select k nearest neighbors (exclude self)
    neighbors <- order(distances)[2:(k + 1)]
    
    # Compute Gaussian weights
    W_gaussian_dense[i, neighbors] <- exp(-distances[neighbors]^2 / (2 * sigma^2))
  }
  
  # Normalize rows to sum to 1 (handle zero rows)
  row_sums <- rowSums(W_gaussian_dense)
  non_zero_rows <- row_sums > epsilon
  W_gaussian_dense[non_zero_rows, ] <- sweep(W_gaussian_dense[non_zero_rows, ], 1, row_sums[non_zero_rows], FUN = "/")
  
  # Calculate sparsity (zero inflation rate)
  zero_count <- sum(W_gaussian_dense == 0)
  total_count <- n_cells * n_cells
  zero_inflation_rate <- zero_count / total_count
  
  # Report sparsity
  message(sprintf("✅ Gaussian weights computed - Zero Inflation Rate (Sparsity): %.2f%%", zero_inflation_rate * 100))
  
  return(W_gaussian_dense)
}

# -------------------------- Core Analysis Function --------------------------
#' Run Complete Baydubig Pipeline (Preprocessing + MCMC + BFDR)
#' 
#' @param raw_expression Raw gene expression matrix (cells × genes)
#' @param raw_coordinates Raw cell coordinates (matrix/data frame, cells × X/Y)
#' @param gene_group_list Named list of gene groups (names = gene names, values = group IDs)
#' @param X Covariate matrix (cells × covariates)
#' @param sigma Gaussian bandwidth for weight calculation (default = 0.01)
#' @param k Number of nearest neighbors for weights (default = 5)
#' @param mcmc_iter Number of MCMC iterations (default = 100)
#' @param mcmc_burn Number of burn-in iterations (default = 50)
#' @param fdr_threshold BFDR control threshold (default = 0.05)
#' @return Integrated result list (preprocessing outputs + MCMC + BFDR results)
run_baydubig <- function(
    raw_expression,
    raw_coordinates,
    gene_group_list,
    X,  # Input raw covariate X (not XWX)
    sigma = 0.01,
    k = 8,
    mcmc_iter = 100,
    mcmc_burn = 50,
    target_bfdr = 0.05,  # BFDR target (0.05 per your definition)
    c_candidates = seq(0.05, 0.9, by = 0.05),  # c values to test
    svg_indices = 1:50  # True SVG indices (for validation only)
) {
  # ===================== Step 1: Preprocessing (Calculate XWX from X + coordinates) =====================
  message("\n===================== Preprocessing =====================")
  # 1.1 Normalize expression
  normalize_expression <- function(expression_matrix) {
    total_counts <- rowSums(expression_matrix)
    mean_expr <- colMeans(expression_matrix)
    var_expr <- apply(expression_matrix, 2, var)
    fit <- lm(var_expr ~ I(mean_expr) + I(mean_expr^2))
    alpha_est <- coef(fit)[3]
    if (alpha_est <= 1e-12) alpha_est <- 1
    anscombe_transformed <- log(expression_matrix + 1 / alpha_est)
    centered_matrix <- sweep(anscombe_transformed, 2, colMeans(anscombe_transformed), FUN = "-")
    return(centered_matrix)
  }
  normalized_expr <- normalize_expression(raw_expression)
  
  # 1.2 Normalize coordinates (Z-score)
  normalize_coordinates <- function(cell_coordinates) {
    cell_coords_mat <- as.matrix(cell_coordinates)
    normalized_coords <- scale(cell_coords_mat, center = TRUE, scale = TRUE)
    rownames(normalized_coords) <- rownames(cell_coordinates)
    colnames(normalized_coords) <- colnames(cell_coordinates)
    return(normalized_coords)
  }
  normalized_coords <- normalize_coordinates(raw_coordinates)
  X_coord <- normalized_coords[, 1]  
  Y_coord <- normalized_coords[, 2]  
  
  # 1.3 Compute Gaussian weights (from coordinates)
  compute_gaussian_weights <- function(cell_coordinates, sigma = 0.01, k = 5) {
    n_cells <- nrow(cell_coordinates)
    W_gaussian_dense <- matrix(0, nrow = n_cells, ncol = n_cells)
    rownames(W_gaussian_dense) <- rownames(cell_coordinates)
    colnames(W_gaussian_dense) <- rownames(cell_coordinates)
    
    for (i in 1:n_cells) {
      coord_i <- cell_coordinates[i, ]
      distances <- sqrt(rowSums((t(t(cell_coordinates) - coord_i))^2))
      neighbors <- order(distances)[2:(k + 1)]  # Exclude self
      W_gaussian_dense[i, neighbors] <- exp(-distances[neighbors]^2 / (2 * sigma^2))
    }
    
    row_sums <- rowSums(W_gaussian_dense)
    non_zero_rows <- row_sums > 1e-2
    W_gaussian_dense[non_zero_rows, ] <- sweep(W_gaussian_dense[non_zero_rows, ], 1, row_sums[non_zero_rows], FUN = "/")
    return(W_gaussian_dense)
  }
  gaussian_weights_dense <- compute_gaussian_weights(normalized_coords, sigma = sigma, k = k)
  weights_sparse <- Matrix::Matrix(gaussian_weights_dense, sparse = TRUE)
  
  # 1.4 Calculate XWX (lag effect) DYNAMICALLY from X + weights
  WX <- gaussian_weights_dense %*% X  # Lag effect of raw covariate X
  XWX <- cbind(X, WX)  # Final covariate matrix (X + lag)
  colnames(XWX) <- c("Covariate_1", "Lag_Covariate_1")
  message("✅ XWX (lag effect) calculated dynamically from X + spatial weights")
  
  # ===================== Step 2: Validate Data =====================
  message("\n===================== Data Validation =====================")
  if (nrow(normalized_expr) != nrow(XWX)) {
    stop("Mismatch in number of cells: expression (", nrow(normalized_expr), ") vs. XWX (", nrow(XWX), ")")
  }
  if (length(gene_group_list) != ncol(normalized_expr)) {
    stop("Mismatch in number of genes: group list (", length(gene_group_list), ") vs. expression (", ncol(normalized_expr), ")")
  }
  
  # ===================== Step 3: Run MCMC Analysis =====================
  message("\n===================== MCMC Analysis =====================")
  # Core MCMC call
  start_time <- Sys.time()
  result_mcmc <- Baydubig_mcmc(
    as.matrix(normalized_expr),
    gene_group_list,
    XWX,  # Use dynamically calculated XWX
    X_coord,
    Y_coord,
    weights_sparse,
    iter = mcmc_iter,
    burn = mcmc_burn
  )
  runtime <- difftime(Sys.time(), start_time, units = "secs")
  message(paste0("✅ MCMC completed - Runtime: ", round(runtime, 2), " seconds"))
  
  # ===================== Step 4: Extract tau_gamma_upper (PPI_j) =====================
  selected_z <- if (!is.null(result_mcmc$selected_Z)) result_mcmc$selected_Z else NA
  tau_gamma_upper <- if (!is.null(result_mcmc[[paste0(selected_z, "_result")]]$tau_gamma_upper)) {
    result_mcmc[[paste0(selected_z, "_result")]]$tau_gamma_upper
  } else {
    stop("tau_gamma_upper (PPI_j) is empty - cannot calculate SVG status")
  }
  if (!all(is.na(tau_gamma_upper))) names(tau_gamma_upper) <- colnames(normalized_expr)
  
  # ===================== Step 5: BFDR Calculation (SVG Status + Gene Names) =====================
  message("\n===================== BFDR Calculation (SVG Status Only) =====================")
  # Load required packages
  if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
  library(dplyr)
  library(tibble)
  
  # Step 5.1: Prepare PPI data (tau_gamma_upper = PPI_j)
  p <- length(tau_gamma_upper)
  ppi_data <- tibble(
    Gene = names(tau_gamma_upper),
    PPI_j = as.numeric(tau_gamma_upper)
  )
  
  # Step 5.2: BFDR(c) calculation (your exact formula)
  calculate_bfdr <- function(ppi_data, c) {
    one_minus_PPI <- 1 - ppi_data$PPI_j
    indicator <- as.numeric(one_minus_PPI < c)
    numerator <- sum(one_minus_PPI * indicator)
    denominator <- sum(indicator)
    bfdr <- ifelse(denominator == 0, Inf, numerator / denominator)
    return(list(bfdr = bfdr, indicator = indicator))
  }
  
  # Step 5.3: Find optimal c (smallest c with BFDR ≤ target)
  bfdr_values <- sapply(c_candidates, function(c) calculate_bfdr(ppi_data, c)$bfdr)
  valid_c_idx <- which(bfdr_values <= target_bfdr)
  
  if (length(valid_c_idx) == 0) {
    stop(paste("No c value found with BFDR ≤", target_bfdr, "- cannot determine SVG status"))
  }
  optimal_c <- c_candidates[min(valid_c_idx)]  # Smallest valid c
  
  # Step 5.4: Determine SVG status (PPI_j ≥ 1 - optimal_c → 1-PPI_j < optimal_c)
  one_minus_PPI <- 1 - ppi_data$PPI_j
  ppi_data <- ppi_data %>%
    mutate(
      is_SVG = one_minus_PPI < optimal_c  # TRUE = SVG, FALSE = non-SVG
    )
  
  # ===================== Step 6: Prepare All Requested Outputs =====================
  # 1. tau_gamma_upper (full PPI results, named vector)
  tau_gamma_results <- ppi_data$PPI_j
  names(tau_gamma_results) <- ppi_data$Gene
  
  # 2. SVG status (named vector: Gene → TRUE/FALSE)
  svg_status <- ppi_data$is_SVG
  names(svg_status) <- ppi_data$Gene
  
  # 3. SVG gene names (only genes with is_SVG = TRUE)
  svg_gene_names <- ppi_data$Gene[ppi_data$is_SVG]
  
  # Final summary
  message("\n🎉 All Results Calculated!")
  message(paste("- Optimal c selected:", optimal_c, "(BFDR ≤", target_bfdr, ")"))
  message(paste("- Total SVGs identified:", length(svg_gene_names), "/", length(svg_status)))
  message("- Output includes: tau_gamma_results, svg_status, svg_gene_names")
  
  # Return structured list with ALL requested outputs
  final_output <- list(
    tau_gamma_results = tau_gamma_results,  # Full tau_gamma/PPI values (named vector)
    svg_status = svg_status,                # Gene → TRUE/FALSE (named vector)
    svg_gene_names = svg_gene_names,        # List of SVG gene names (character vector)
    optimal_c = optimal_c,                  # Selected c (for reference)
    target_bfdr = target_bfdr               # Target BFDR (for reference)
  )
  
  return(final_output)
}
