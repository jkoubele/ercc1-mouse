library(Seurat)
library(tidyverse)


options(future.globals.maxSize = 4 * 1024^3)  # Increases parallelization limit (used by SCTransform) to 4 GB

path_prefix <- "/ercc1/data"

input_folder <- file.path(path_prefix, "cell_quality_control_output")
output_folder <- file.path(path_prefix, "clustering_output")
dir.create(output_folder)

seurat_object <- readRDS(file.path(input_folder, "seurat_object.rds")) |>
  SCTransform() |>
  RunPCA()

# PCA diagnostics plots
pca_by_sample <- DimPlot(seurat_object, reduction = 'pca', group.by = 'sample') +
  labs(title = "PCA by sample")
ggsave(file.path(output_folder, "pca_sample.png"),
       plot = pca_by_sample,
       bg = 'white')

elbow_plot <- ElbowPlot(seurat_object, ndims = 50) +
  labs(title = "PCA Elbow Plot")
ggsave(file.path(output_folder, "pca_elbow_plot.png"),
       plot = elbow_plot,
       bg = 'white')

# Run clustering and UMAP
num_pca_dims <- 30 # based on the SCTransform vignette https://satijalab.org/seurat/articles/sctransform_vignette.html
seurat_object <- seurat_object |>
  FindNeighbors(dims = 1:num_pca_dims, reduction = 'pca') |>
  FindClusters() |>
  RunUMAP(dims = 1:num_pca_dims)

# UMAP plots
umap_by_sample <- DimPlot(seurat_object, reduction = 'umap', group.by = "sample") +
  labs(title = "UMAP by sample")
ggsave(file.path(output_folder, "umap_sample.png"),
       plot = umap_by_sample,
       bg = 'white')

umap_by_cluster <- DimPlot(seurat_object, reduction = 'umap', group.by = "seurat_clusters") +
  labs(title = "UMAP by cluster")
ggsave(file.path(output_folder, "umap_cluster.png"),
       plot = umap_by_cluster,
       bg = 'white')

# PCA plot grouped by clusters
pca_by_sample <- DimPlot(seurat_object, reduction = 'pca', group.by = 'seurat_clusters') +
  labs(title = "PCA by clusters")
ggsave(file.path(output_folder, "pca_cluster.png"),
       plot = pca_by_sample,
       bg = 'white')

saveRDS(seurat_object, file.path(output_folder, 'seurat_object.rds'))

