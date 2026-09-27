# BayDuBiG

**Bayesian Durbin-based Bi-level Identification of Spatially Variable Genes**

BayDuBiG is a Bayesian bi-level variable-selection model for identifying
**spatially variable genes (SVGs)** in spatial transcriptomics data. It couples a
spatial autoregressive (Durbin) term with a two-level spike-and-slab prior,
jointly modeling spatial spillover effects and coordinated regulation within
functional gene groups.

This repository accompanies the paper *Bayesian Durbin-based bi-level
identification of spatially variable genes* .

## Model overview

- **Spatial autoregressive (Durbin) term**: a spatial weight matrix $\mathbf{W}$
  and a gene-specific autoregressive coefficient $\rho_j$ capture spatial
  dependence and spillover among neighboring spots; covariate-induced spillover
  is modeled through the spatially lagged covariates $\mathbf{W}\mathbf{X}$.
- **Two-level spike-and-slab prior**:
  - **group (pathway) level** $\tau$: whether a pathway is associated with spatial variation;
  - **within-group (gene) level** $\gamma$: whether a gene inside a pathway is an SVG.
- **Union semantics for the group effect**: a gene is flagged as significant if
  it belongs to **at least one** selected group (pathway):

  $$\zeta_j = 1 - \prod_{g \in \mathcal{G}_j}(1 - \tau_g)$$

  where $\mathcal{G}_j$ is the set of groups gene $j$ belongs to.
- **BFDR control**: a Bayesian false discovery rate (BFDR) rule is applied to the
  posterior inclusion probabilities (PPIs) to control the overall error rate
  (default target: 0.05).

## Repository contents

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
└── demo/                  # a ready-to-run example on a real-data slice
    ├── demo_data.rds          # demo data: a real MBSA slice (400 spots x 1000 genes, 339 pathways)
    ├── demo_expr_4genes.png   # figure: expression of 4 genes across the slice
    ├── demo_covariate_pie.png # figure: cell-type composition at each spot
    └── Demo.R                 # demo script (loads data, draws figures, runs the pipeline)
```

- **`README.md`** — the starting point.
- **`BayDuBiG/`** — the R package that implements the method.
- **`demo/`** — a self-contained example to test the method right away.

## Installation

Required R packages: `Rcpp`, `RcppArmadillo`, `Matrix`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "Matrix"))
```

**Option 1 — install as an R package** (from the repository root):

```bash
R CMD INSTALL BayDuBiG
```

```r
library(BayDuBiG)
```

**Option 2 — use without installing** (from the repository root):

```r
source("BayDuBiG/R/BayDuBiG.R")
Rcpp::sourceCpp("BayDuBiG/src/BayDuBiGMCMC.cpp")
```

## Usage

```r
results <- run_BayDuBiG(
  raw_expression    = expression,   # gene expression (cells x genes)
  raw_coordinates   = coordinates,  # spatial coordinates (cells x 2)
  gene_group_list   = groups,       # gene -> pathway index(es), starting from 1
  X                 = X,            # spot covariates (cells x q), or NULL
  iter              = 200,
  burn              = 100,
  target_bfdr       = 0.05,
  informative_gamma = TRUE          # informative prior used in the paper
)

results$svg_gene_names   # names of the identified SVGs
results$gene_results     # per-gene table: PPI, basis log-likelihoods, selected basis, SVG flag
```

## Running the demo

A ready-to-run example on a real MBSA data slice is in `demo/`. From the
repository root:

```r
source("demo/Demo.R")
```

It reads `demo/demo_data.rds` (400 spots x 1000 genes with real KEGG pathway
annotations), draws two diagnostic figures, and runs the pipeline. Gene and spot
names are anonymized to `Gene_1..Gene_1000` / `Spot_1..Spot_400`.

## OpenMP acceleration

The MCMC core supports OpenMP for the group-level $\tau$ updates. It is disabled
by default (Apple clang lacks `-fopenmp`); to enable it, point
`BayDuBiG/src/Makevars` (or `~/.R/Makevars`) to an OpenMP-capable compiler, then
set the number of threads before running the model:

```r
Sys.setenv(OMP_NUM_THREADS = 8)   # adjust to your CPU cores
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
