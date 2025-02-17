library(Seurat)
library(tidyverse)
library(rjson)


path_prefix <- "/ercc1/data"

input_folder <- file.path(path_prefix, "scDblFinder_output")
output_folder <- file.path(path_prefix, "cell_quality_control_output")
dir.create(output_folder)

seurat_object <- readRDS(file.path(input_folder, "seurat_object.rds")) |>
  PercentageFeatureSet("^mt-", col.name = "percent_mito")

scatter_count_vs_genes <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme_minimal()
ggsave(file.path(output_folder, "scatter_count_vs_genes.png"), plot = scatter_count_vs_genes, bg = 'white')

qc_by_metric <- function(seurat_object, metric_name, num_mads = 3) {
  # Adapted from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-cells
  # For comparison, please note that scipy.stats.median_abs_deviation doesn't scale the result, while
  # the R function mad() scales the output by a constant of 1.4826
  metric <- seurat_object@meta.data[[metric_name]]
  lower_threshold <- median(metric) - num_mads * mad(metric)
  upper_threshold <- median(metric) + num_mads * mad(metric)
  return(
    list(
      is_outlier = (metric < lower_threshold) | (metric > upper_threshold),
      lower_threshold = lower_threshold,
      upper_threshold = upper_threshold
    )
  )
}


seurat_object@meta.data$keep_after_qc <- seurat_object@meta.data$scDblFinder.class == 'singlet'

for (metric_name in c('nCount_RNA', 'nFeature_RNA')) {
  print(metric_name)
  qc_info <- qc_by_metric(seurat_object, metric_name)
  histogram <- ggplot(seurat_object@meta.data,
                      aes(x = .data[[metric_name]])) +
    geom_histogram() +
    theme_minimal() +
    labs(title = paste("Histogram of", metric_name), x = metric_name, y = "Count") +
    geom_vline(xintercept = qc_info$upper_threshold, color = "red")
  if (qc_info$lower_threshold > 0) {
    histogram <- histogram + geom_vline(xintercept = qc_info$lower_threshold, color = "red")
  }

  print(histogram)
  ggsave(file.path(output_folder, paste0("histogram_of_", metric_name, ".png")),
         plot = histogram,
         bg = 'white')

  seurat_object@meta.data$keep_after_qc <- seurat_object@meta.data$keep_after_qc * (qc_info$is_outlier == FALSE)
}

# Over 91% of cells have 0 mitochondrial genes detected, we thus set threshold manually to low but permissive value
print(paste("Fraction of cells with 0 mitochondiral genes:", mean(seurat_object@meta.data$percent_mito == 0)))
seurat_object@meta.data$keep_after_qc <- seurat_object@meta.data$keep_after_qc * (seurat_object@meta.data$percent_mito < 0.1)

write(
  toJSON(
    list(
      cells_before_qc = length(seurat_object@meta.data$keep_after_qc),
      cells_after_qc = sum(seurat_object@meta.data$keep_after_qc),
      fraction_of_cells_filtered = 1 - mean(seurat_object@meta.data$keep_after_qc),
      fraction_of_doublets = mean(seurat_object@meta.data$scDblFinder.class == 'doublet')
    )
  ),
  file = file.path(output_folder, "qc_filtering_stats.json"))

seurat_object <- subset(seurat_object, subset = keep_after_qc == TRUE)
saveRDS(seurat_object, file.path(output_folder, 'seurat_object.rds'))


