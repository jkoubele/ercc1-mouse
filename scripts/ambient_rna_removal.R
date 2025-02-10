library(Seurat)
library(SingleCellExperiment)
library(decontX)
library(rjson)


path_prefix <- "/ercc1/data"
input_folder <- file.path(path_prefix, "cellranger_output")
output_folder <- file.path(path_prefix, "decontX_output")
dir.create(output_folder)

fraction_of_umis_removed = list()
for (sample_folder in list.dirs(input_folder, recursive = FALSE)) {
  sample_name <- basename(sample_folder)
  print(sample_name)
  expression_matrix_filtered <- Read10X(
    file.path(sample_folder, "outs/filtered_feature_bc_matrix/"))

  expression_matrix_raw <- Read10X(file.path(sample_folder, "outs/raw_feature_bc_matrix/"))

  sce_filtered <- SingleCellExperiment(list(counts = expression_matrix_filtered))
  sce_raw <- SingleCellExperiment(list(counts = expression_matrix_raw))
  sce_decont <- decontX(sce_filtered, background = sce_raw)

  seurat_object_decont <- CreateSeuratObject(round(decontXcounts(sce_decont)))

  fraction_of_umis_removed[[sample_name]] <- 1 - sum(seurat_object_decont[['RNA']]$counts) / sum(expression_matrix_filtered)

  output_subfolder <- file.path(output_folder, sample_name)
  saveRDS(seurat_object_decont, file.path(output_subfolder, 'seurat_object.rds'))

}

write(toJSON(fraction_of_umis_removed), file = file.path(output_folder, "fraction_of_umis_removed.json"))

