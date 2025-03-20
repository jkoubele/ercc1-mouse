library(Seurat)
library(rjson)
library(tidyverse)
library(Polychrome)


source("https://raw.githubusercontent.com/jkoubele/sc-type-refactored/main/cell_type_identification.R")

path_prefix <- "/cellfile/datapublic/jkoubele/ercc1/data/"
# path_prefix <- "/ercc1/data/"
input_folder <- file.path(path_prefix, "clustering_output")
output_folder <- file.path(path_prefix, "cell_type_annotation")

seurat_object <- readRDS(file.path(input_folder, "seurat_object.rds"))

mouse_markers <- fromJSON(file = "https://raw.githubusercontent.com/jkoubele/sc-type-refactored/main/cell_type_markers_mouse.json")

kidney_markers <- mouse_markers$Kidney

cell_type_scores <- cell_type_scoring(seurat_object, kidney_markers)
clusters_cell_type <- cluster_cell_type_clasification(seurat_object, cell_type_scores)
seurat_object$cell_type <- clusters_cell_type$cell_type[match(seurat_object$seurat_clusters, clusters_cell_type$seurat_clusters)]

umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "cell_type", shuffle = TRUE) +
  labs(title = "Cell type UMAP", color = "Cell type") +  # Use `color =` for discrete groups
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 15)
  )

print(umap_plot)
ggsave(file.path(output_folder, "cell_type_umap.png"),
       plot = umap_plot,
       bg = 'white')

# Cell type counts 
cell_type_count_barplot <- ggplot(seurat_object@meta.data, aes(y = fct_rev(fct_infreq(cell_type)))) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(y = "Cell Type",
       x = "Number of Cells",
       title = "Cell type counts") +
  theme(axis.text.x = element_text(12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 20, hjust = 0.5))

print(cell_type_count_barplot)
ggsave(file.path(output_folder, "cell_type_count_barplot.png"),
       plot = cell_type_count_barplot,
       bg = 'white')

# Sample composition plot 
sample_composition_data <- seurat_object@meta.data |>
  count(sample_name, cell_type) |>
  group_by(sample_name) |>
  mutate(proportion = n / sum(n))

polychrome_seed <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3", "#00CED1")
color_palette <- createPalette(length(sample_composition_data$cell_type), seedcolors = polychrome_seed)
names(color_palette) <- unique(sample_composition_data$cell_type)

sample_composition_plot <- ggplot(sample_composition_data, aes(x = sample_name, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 11)) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Sample",
       y = "Proportion",
       fill = "Cell Type",
       title = "Samples cell type composition")

print(sample_composition_plot)

ggsave(file.path(output_folder, "sample_composition_plot.png"),
       plot = sample_composition_plot,
       bg = 'white')

saveRDS(seurat_object, file.path(output_folder, 'seurat_object.rds'))
