# ==============================================================================
# Run Baydubig Pipeline on demo_data.rds (20x20 Grid, 100 Genes Example)
# ==============================================================================


source(here("Baydubig","R", "BayDuBiG.R")) 
compile_baydubig_cpp(cpp_path = here("Baydubig", "src", "Baydubig_mcmc.cpp"))




# -------------------------- Example Usage --------------------------
demo_data <- readRDS(here("demo_data.rds"))

results <- run_baydubig(
  raw_expression = demo_data$raw_expression,
  raw_coordinates = demo_data$raw_coordinates,
  gene_group_list = demo_data$gene_group_list,
  X = demo_data$X,
  sigma = demo_data$simulation_params$sigma,
  k = demo_data$simulation_params$k,
  mcmc_iter = 1000,
  mcmc_burn = 500,
  target_bfdr = 0.05,
  c_candidates = seq(0.05, 0.9, by = 0.05),
  svg_indices = demo_data$svg_indices
)

svg_taugamma_df <- data.frame(
  Gene = names(results$tau_gamma_results),
  tau_gamma = results$tau_gamma_results,
  is_SVG = results$svg_status,
  row.names = NULL
)
