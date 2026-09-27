# BayDuBiG

**Bayesian Durbin-based Bi-level Identification of Spatially Variable Genes**

BayDuBiG is a Bayesian bi-level variable-selection model for identifying
**spatially variable genes (SVGs)** in spatial transcriptomics data. It couples a
spatial autoregressive (Durbin) term with a two-level spike-and-slab prior, so
that both spatial spillover effects and coordinated regulation within functional
gene groups are explicitly modeled.

This repository accompanies the paper *Bayesian Durbin-based bi-level
identification of spatially variable genes* (Feng, Wang, and Wu).

## Model overview

- **Spatial autoregressive (Durbin) term**: a spatial weight matrix `W` and a
  gene-specific autoregressive coefficient `rho_j` capture spatial dependence and
  spillover among neighboring spots. Covariate-induced spillover is modeled
  through the spatially lagged covariates `W X`.
- **Two-level spike-and-slab prior**:
  - **group (pathway) level** `tau`: whether a pathway is associated with spatial variation;
  - **within-group (gene) level** `gamma`: whether a gene inside a pathway is an SVG.
- **Union semantics for the group effect**: a gene is flagged as significant if
  it belongs to **at least one** selected group (pathway):

  $$\zeta_j = 1 - \prod_{g \in \mathcal{G}_j}(1 - \tau_g)$$

  where $\mathcal{G}_j$ is the set of groups gene $j$ belongs to.
- **BFDR control**: a Bayesian false discovery rate (BFDR) rule is applied to the
  posterior inclusion probabilities (PPIs) to control the overall error rate
  (default target: 0.05).

## Repository contents

When you open this repository on GitHub you will see three parts:

```
.
├── README.md              # this file (the entry point)
├── .gitignore
├── BayDuBiG/              # the R package: source code of the method
│   ├── DESCRIPTION        #   package metadata
│   ├── NAMESPACE          #   namespace
│   ├── R/
│   │   ├── BayDuBiG.R     #     R entry points (preprocessing + run_BayDuBiG + BFDR)
│   │   └── RcppExports.R  #     auto-generated C++ interface
│   └── src/
│       ├── BayDuBiGMCMC.cpp  #  C++ MCMC core (RcppArmadillo, optional OpenMP)
│       ├── RcppExports.cpp   #  auto-generated C++ registration
│       └── Makevars          #  build configuration
└── demo/                  # a self-contained, ready-to-run example
    ├── demo_data.rds          # demo data: a real MBSA slice (400 spots x 1000 genes, 339 pathways)
    ├── demo_expr_4genes.png   # figure: expression of 4 genes across the slice
    ├── demo_covariate_pie.png # figure: cell-type composition at each spot
    └── Demo.R                 # demo script (loads data, draws figures, runs the pipeline)
```

- **`README.md`** — what you are reading; the starting point.
- **`BayDuBiG/`** — the R package that implements the method (installable with
  `R CMD INSTALL`).
- **`demo/`** — a ready-to-run example on a real-data slice, so you can test the
  method immediately without preparing your own data.

## Installation

Required R packages: `Rcpp`, `RcppArmadillo`, `Matrix`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "Matrix"))
```

### Option 1 — install as an R package (recommended)

From a terminal, in the repository root (the directory containing this README):

```bash
R CMD INSTALL BayDuBiG
```

Or from within R:

```r
devtools::install("BayDuBiG")
```

Then load and use it:

```r
library(BayDuBiG)
res <- run_BayDuBiG(...)
```

### Option 2 — use without installing (quick try-out)

From the repository root:

```r
source("BayDuBiG/R/BayDuBiG.R")
Rcpp::sourceCpp("BayDuBiG/src/BayDuBiGMCMC.cpp")
```

## Quick start

```r
source("BayDuBiG/R/BayDuBiG.R")
Rcpp::sourceCpp("BayDuBiG/src/BayDuBiGMCMC.cpp")

# Load your data (convention: rows = cells, columns = genes)
#   expression  : gene expression matrix (cells x genes)
#   coordinates : spatial coordinates (cells x 2)
#   groups      : named list mapping each gene to its pathway index(es)
#                 (indices start at 1; use integer(0) if a gene has no pathway)
#   X           : covariate matrix (cells x covariates), or NULL

results <- run_BayDuBiG(
  raw_expression  = expression,
  raw_coordinates = coordinates,
  gene_group_list = groups,
  X               = X,
  iter            = 2000,
  burn            = 1000,
  target_bfdr     = 0.05
)

results$svg_gene_names    # names of the identified SVGs
results$tau_gamma_results # per-gene PPI
results$svg_status        # per-gene logical flag
results$gene_results      # per-gene table: PPI, the three basis log-likelihoods,
                          # the selected basis, and the SVG flag
results$basis_per_gene    # the basis selected for each gene
```

## Running the demo

The demo is the easiest way to test the method with a self-contained example.
From the repository root, run:

```r
source("demo/Demo.R")
```

This reads `demo/demo_data.rds`, a self-contained slice of the real MBSA spatial
transcriptomics dataset: 400 spots (taken from the bottom-left corner of the
tissue) x 1000 genes with real KEGG pathway annotations (339 pathways), plus the
real cell-type proportions as the covariate `X`. Gene and spot names are
anonymized to `Gene_1..Gene_1000` / `Spot_1..Spot_400`.

Before running the pipeline, the demo draws two diagnostic figures:

- `demo/demo_expr_4genes.png` — expression of the 4 most spatially varying
  genes across the slice;
- `demo/demo_covariate_pie.png` — cell-type composition of the covariate `X`
  at each spot (one pie per spot).

It then runs the full pipeline (with `informative_gamma = TRUE`) and prints the
identified SVGs, the per-gene result table (PPI, three basis log-likelihoods,
selected basis, SVG flag), and the pathway names. To keep the demo fast it uses a
small number of MCMC iterations (`iter = 200`); for a full analysis use
`iter = 2000, burn = 1000`. The figures require the `ggplot2` and `scatterpie`
packages.

## Main parameters

| Parameter | Meaning | Default |
|-----------|---------|---------|
| `sigma` | Gaussian spatial-weight bandwidth | 0.01 |
| `k` | number of nearest neighbors for the Gaussian weights | 8 |
| `iter` | total MCMC iterations | 2000 |
| `burn` | burn-in iterations | 1000 |
| `target_bfdr` | BFDR control target | 0.05 |
| `informative_gamma` | use the informative gamma prior | FALSE |

### Informative gamma prior

When `informative_gamma = TRUE`, the gene-level prior switches to a stronger,
pathway-informed specification: `Beta(10, 1)` when the group effect is active and
`Beta(1, 10)` when it is inactive. This is the setting used in the accompanying
paper; pass `informative_gamma = TRUE` to reproduce it.

## OpenMP acceleration

The MCMC main loop in `BayDuBiG/src/BayDuBiGMCMC.cpp` supports OpenMP. By default
it is compiled with the system compiler (which may be Apple clang, lacking
`-fopenmp`), in which case the code automatically falls back to a serial
implementation with identical results. To enable OpenMP, edit
`BayDuBiG/src/Makevars` (see the comments there) or set `~/.R/Makevars` to point
to an OpenMP-capable compiler before installing.

Once OpenMP is enabled, increasing the number of threads speeds up the
group-level `tau` updates (the main computational bottleneck). Set the number of
threads before running the model, for example:

```r
Sys.setenv(OMP_NUM_THREADS = 8)   # use 8 threads (adjust to your CPU cores)
```

## Citation

If this package is useful for your research, please cite the corresponding paper:

> Feng, X., Wang, T., and Wu, M. *Bayesian Durbin-based bi-level identification
> of spatially variable genes.* (add the journal, volume, and DOI when available)

and the software:

```bibtex
@software{BayDuBiG,
  title   = {BayDuBiG: Bayesian Durbin-based Bi-level Identification of Spatially Variable Genes},
  author  = {Xingdong Feng and Tianyi Wang and Mengyun Wu},
  url     = {https://github.com/mengyunwu2020/BayDuBiG}
}
```

## Authors

Xingdong Feng, Tianyi Wang, and Mengyun Wu — School of Statistics and Data
Science, Shanghai University of Finance and Economics. Correspondence:
Mengyun Wu (wu.mengyun@mail.shufe.edu.cn).
