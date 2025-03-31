library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)


path_prefix <- "/ercc1/data"

input_folder <- file.path(path_prefix, "decontX_output")
output_folder <- file.path(path_prefix, "scDblFinder_output")
dir.create(output_folder)

sample_list <- list()
for (sample_folder in list.dirs(input_folder, recursive = FALSE)) {
  sample_name <- basename(sample_folder)
  print(sample_name)
  sample_seurat_object <- readRDS(file.path(sample_folder, "seurat_object.rds"))
  sample_seurat_object$sample <- sample_name
  sample_seurat_object <- RenameCells(sample_seurat_object, add.cell.id = sample_name)
  sample_list[[sample_name]] <- sample_seurat_object
}

merged_seurat_object <- merge(sample_list[[1]],
                              y = sample_list[-1]) |>
  JoinLayers() |>
  subset(subset = nCount_RNA > 200) # removal of ambient RNA may left some cells with very low UMI count
# we apply very permissive filtering for now, as recommended by scDblFinder warning

sce <- as.SingleCellExperiment(merged_seurat_object) |>
  scDblFinder(samples = "sample",
              verbose = TRUE,
              BPPARAM = MulticoreParam(4))

seurat_with_detected_doublets <- as.Seurat(sce, data = NULL)
saveRDS(seurat_with_detected_doublets, file.path(output_folder, 'seurat_object.rds'))
