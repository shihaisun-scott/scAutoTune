# sample script to analyse the PBMC dataset

# libraries
library(scAutoTune)

# load the file
obj <- readRDS("../inst/sample_data/pbmc_cluster0.rds")

# analyse using scAutoTune
# 1. find the optimum number of pcs
#feature_steps <- seq(1000, 3000, by = 1000)
feature_steps = seq(500, 6000, by = 500)
max_pcs <- 20
pc_out <- autotune_find_pcs(obj,
                            feature_steps = feature_steps,
                            max_pcs = max_pcs)
pc_out$suggested_pcs
pc_out$plot

# 2. find the optimum number of features and resolution
res <- autotune_find_features_resolution(
  obj,
  n_pcs = 5,
  nfeatures_range =  seq(500, 3000, by = 50),
  resolutions = seq(0.2, 1.0, by = 0.02),
  k_val = 10,
  output_csv = "cluster0_sweep.csv",
  output_pdf = "cluster0_heatmap_umap.pdf"
)

res$chosen_params
print(res$heatmap)
print(res$umap_cluster_plot)
