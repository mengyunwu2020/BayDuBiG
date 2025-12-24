# BayDuBiG: Bayesian Dubin-based Bi-level framework for identifying SV Genes

BayDuBiG is a **Bayesian hierarchical model** (inspired by the Durbin model) designed to identify **Spatially Variable Genes (SVGs)** from spatial transcriptomics data. It uses a **bi-level (group-gene) inference framework** combined with MCMC sampling (implemented in C++ for efficiency) to detect genes with non-random spatial expression patterns.


## Repository Structure
```
./
├── BayDubig/                # Core code main directory
│   ├── R/                   # R wrapper code (preprocessing, pipeline call, result parsing)
│   │   ├── BayDuBiG.R       # Main pipeline function (run_baydubig)
│   │   └── RcppExports.R    # Auto-generated Rcpp binding file (do not modify)
│   ├── src/                 # C++ core code (MCMC sampling implementation)
│   │   ├── Baydubig_mcmc.cpp # Core MCMC sampling logic (OpenMP parallel acceleration)
│   │   └── RcppExports.cpp  # Auto-generated Rcpp binding file (do not modify)
│   ├── NAMESPACE            # R package namespace definition (exports public functions)
│   └── DESCRIPTION          # R package metadata (version, dependencies, author, etc.)
├── Demo/                    # Demo directory (quick start)
│   ├── Demo.R               # Example running script (complete call workflow)
│   ├── demo_data.rds        # Example dataset (100 genes × 400 cells, 20×20 spatial grid)
└── README.md                # Repository main documentation (this file)
```


## Core Functionality
BayDuBiG implements:
1. **Preprocessing**: Expression normalization (Anscombe transform) + coordinate Z-score standardization + Gaussian spatial weight calculation.
2. **Bi-level Bayesian Inference**:
   - Group-level: Shared spatial patterns across gene groups.
   - Gene-level: Individual gene spatial variability.
3. **MCMC Sampling**: Efficient C++ implementation (OpenMP-parallelized) for posterior inference.
4. **BFDR Control**: False discovery rate correction to select significant SVGs.


## Installation
### Prerequisites
- R with packages: `Rcpp`, `RcppEigen`, `Matrix`, `dplyr`, `tibble`, `here`
- C++ compiler 


## Quick Start (Using the Demo)
The `Demo/Demo.R` script runs the full pipeline on the example dataset (`demo_data.rds`):

1. Navigate to the `Demo/` folder.
2. Run `Demo.R` in R:
   ```r
   # Load required packages
   library(here)
   source(here("..", "R", "BayDuBiG.R"))
   
   # Load example data (100 genes × 400 cells, 20x20 grid)
   demo_data <- readRDS(here("demo_data.rds"))
   
   # Run BayDuBiG pipeline
   results <- run_baydubig(
     raw_expression = demo_data$raw_expression,
     raw_coordinates = demo_data$raw_coordinates,
     gene_group_list = demo_data$gene_group_list,
     X = demo_data$X,          # Raw covariate matrix
     sigma = 0.01,             # Gaussian weight bandwidth
     k = 8,                    # Number of nearest neighbors
     mcmc_iter = 100,          # MCMC iterations
     mcmc_burn = 50,           # Burn-in iterations
     target_bfdr = 0.05        # BFDR threshold for SVG selection
   )
   
   # Access results
   print("tau_gamma (PPI) values for genes:")
   print(head(results$tau_gamma_results))
   
   print("SVG status (TRUE = SVG):")
   print(head(results$svg_status))
   
   print("Identified SVG gene names:")
   print(head(results$svg_gene_names))
   ```


## Key Outputs
The `run_baydubig` function returns a structured list with:
- `tau_gamma_results`: Named vector of PPI (Posterior Probability of Importance) values for each gene.
- `svg_status`: Named logical vector (TRUE = gene is an SVG).
- `svg_gene_names`: Character vector of genes identified as SVGs.
- `optimal_c`: Threshold `c` used for BFDR control.


## Method Overview
BayDuBiG uses a **bi-level hierarchical model** to model spatial gene expression:
1. **Spatial Structure**: Encodes spatial dependence via a Gaussian weight matrix (k-nearest neighbors).
2. **Group-Gene Hierarchy**: Shares information across gene groups to improve power for rare SVGs.
3. **MCMC Inference**: Samples posterior distributions of spatial parameters (ρ, τ, γ) to quantify gene-level spatial variability.
4. **BFDR Control**: Selects SVGs by controlling the Bayesian False Discovery Rate (BFDR ≤ 0.05 by default).


