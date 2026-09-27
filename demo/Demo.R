# ==============================================================================
# BayDuBiG demo — a real MBSA data slice
#
# Loads demo/demo_data.rds (a self-contained slice of the real MBSA dataset),
# draws two quick diagnostic figures, then runs the full BayDuBiG pipeline.
#
# To run (working directory = repository root):
#   source("BayDuBiG/R/BayDuBiG.R")
#   Rcpp::sourceCpp("BayDuBiG/src/BayDuBiGMCMC.cpp")
#   source("demo/Demo.R")
#
# Requires (for the figures): ggplot2, scatterpie
# ==============================================================================

source("BayDuBiG/R/BayDuBiG.R")
Rcpp::sourceCpp("BayDuBiG/src/BayDuBiGMCMC.cpp", verbose = FALSE)

# ---------- 1. Load the demo data ----------
demo_data <- readRDS("demo/demo_data.rds")
expr <- demo_data$raw_expression
coords <- demo_data$raw_coordinates
X <- demo_data$X

cat("\n==================== Demo data ====================\n")
cat("Spots:", nrow(expr),
    " | Genes:", ncol(expr),
    " | Groups:", length(unique(unlist(demo_data$gene_group_list))), "\n")
cat("Covariate X (cell-type proportions):", paste(dim(X), collapse = "x"), "\n")

# ---------- 2. Figure 1: expression of 4 genes on the slice ----------
suppressMessages(library(ggplot2))
vars <- apply(expr, 2, stats::var)
top4 <- names(sort(vars, decreasing = TRUE)[1:4])   # 4 most spatially varying genes

plot_df <- do.call(rbind, lapply(top4, function(g) {
  data.frame(x = coords[, 1], y = coords[, 2], expression = expr[, g], gene = g)
}))
p1 <- ggplot(plot_df, aes(x = x, y = y, color = expression)) +
  geom_point(size = 1.2) +
  scale_color_gradient(low = "#e5f5e0", high = "#006d2c") +
  facet_wrap(~ gene, nrow = 1) +
  coord_fixed() +
  labs(title = "Expression of 4 genes across the tissue slice") +
  theme_minimal()
print(p1)

# ---------- 3. Figure 2: covariate (cell-type) composition pie per spot ----------
suppressMessages(library(scatterpie))
ct_short <- c("Type.I.SGN", "Neuron", "Neuroendocrine", "Astrocyte",
              "Oligodendrocyte", "Goblet", "Tuft")   # short labels for the legend
pie_df <- data.frame(x = coords[, 1], y = coords[, 2], as.data.frame(X),
                     check.names = FALSE)
colnames(pie_df)[3:9] <- ct_short
p2 <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = 80),
                  data = pie_df, cols = ct_short, color = NA) +
  coord_fixed() +
  labs(title = "Cell-type composition at each spot") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 11))
print(p2)

# ---------- 4. Run the full BayDuBiG pipeline ----------
set.seed(1)
results <- run_BayDuBiG(
  raw_expression    = expr,
  raw_coordinates   = coords,
  gene_group_list   = demo_data$gene_group_list,
  X                 = X,
  iter              = 200,
  burn              = 100,
  target_bfdr       = 0.05,
  informative_gamma = TRUE
)

cat("\n==================== Results ====================\n")
cat("Covariate basis selected per gene:\n")
print(table(results$basis_per_gene))
cat("Number of SVGs identified:", length(results$svg_gene_names),
    "/", ncol(expr), "\n")

cat("\nPer-gene results (top 10 by PPI):\n")
gr <- results$gene_results
gr <- gr[order(-gr$tau_gamma), ]
print(head(gr, 10), row.names = FALSE)

cat("\nGroup information (first 5 groups):\n")
print(head(demo_data$group_info, 5), row.names = FALSE)

cat("\nThe full per-gene table is available as results$gene_results.\n")
