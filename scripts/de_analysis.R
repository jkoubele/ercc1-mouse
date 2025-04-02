library(Seurat)
library(tidyverse)
library(topGO)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


path_prefix <- "/cellfile/datapublic/jkoubele/ercc1/data/"
# path_prefix <- "/ercc1/data/"

input_folder <- file.path(path_prefix, "cell_type_annotation")
output_folder <- file.path(path_prefix, "de_analysis")
dir.create(output_folder)

seurat_object <- readRDS(file.path(input_folder, "seurat_object.rds"))

gene_name_to_entrez_id <- bitr(rownames(seurat_object),
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Mm.eg.db) |>
  rename(gene = SYMBOL)

save_de_and_go_enrichment <- function(de,
                                      expressed_gene_names,
                                      output_subfolder) {
  expressed_genes <- tibble(gene = expressed_gene_names) |>
    left_join(gene_name_to_entrez_id, by = "gene") |>
    drop_na()

  de_up <- de |>
    dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05)
  de_down <- de |>
    dplyr::filter(avg_log2FC < 0, p_val_adj < 0.05)

  write_tsv(de, file.path(output_subfolder, "de_all.tsv"))
  write_tsv(de_up, file.path(output_subfolder, "de_significant_up.tsv"))
  write_tsv(de_down, file.path(output_subfolder, "de_significant_down.tsv"))


  de_list <- list(up = de_up, down = de_down)
  go_output_folder <- file.path(output_subfolder, 'GO_enrichment')
  dir.create(go_output_folder)

  for (direction in names(de_list)) {
    de_dataframe <- de_list[[direction]] |> drop_na()

    for (go_ontology in c('MF', 'BP', 'CC')) {
      go_enrichment <- enrichGO(gene = de_dataframe$ENTREZID,
                                OrgDb = org.Mm.eg.db,
                                ont = go_ontology,
                                universe = expressed_genes$ENTREZID,
                                readable = TRUE)
      write_tsv(go_enrichment@result, file.path(go_output_folder, paste0('go_', go_ontology, '_', direction, '_all.tsv')))

      go_enrichment_simplified <- go_enrichment |> clusterProfiler::simplify()
      write_tsv(go_enrichment_simplified@result, file.path(go_output_folder, paste0('go_', go_ontology, '_', direction, '_simplified.tsv')))

      top_go_terms <- go_enrichment_simplified@result |>
        arrange(p.adjust) |>
        head(15) |>
        mutate(
          gene_fraction = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
          Description = ifelse(nchar(Description) > 50,
                               paste0(str_sub(Description, 1, 50), '...'),
                               Description)
        ) |>
        mutate(Description = factor(Description, levels = rev(Description)))

      go_dotplot <- ggplot(top_go_terms, aes(x = gene_fraction, y = Description, size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme_minimal() +
        labs(
          title = paste0('GO ', go_ontology, ' Enrichment (', direction, 'regulated genes)'),
          x = "Gene Ratio",
          y = "GO Term",
          color = "Adjusted p-value",
          size = "Num. genes"
        ) +
        scale_size_continuous(range = c(4, 9)) +
        theme(axis.text = element_text(size = 12))

      print(go_dotplot)
      ggsave(file.path(go_output_folder, paste0('GO_', go_ontology, '_', direction, '.png')),
             plot = go_dotplot,
             bg = 'white')


    }

  }

}


# Compare KD vs. Wt using all cells
output_subfolder_all_cells <- file.path(output_folder, 'KD_vs_Wt_all_cells')
dir.create(output_subfolder_all_cells)

gene_sct_count_total <- rowSums(seurat_object@assays$SCT$counts)
expressed_gene_names_all_cells <- names(gene_sct_count_total[gene_sct_count_total > 0])

de <- FindMarkers(seurat_object,
                  ident.1 = 'KD',
                  ident.2 = 'Wt',
                  group.by = 'genotype') |>
  rownames_to_column('gene') |>
  left_join(gene_name_to_entrez_id, by = 'gene')

save_de_and_go_enrichment(de, expressed_gene_names_all_cells, output_subfolder_all_cells)


Idents(seurat_object) <- "cell_type"
for (cell_type in unique(seurat_object$cell_type)) {

  if ((!'Wt' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type])) |
    (!'KD' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type]))) {
    next # cell type is not present in both genotype groups
  }


  de <- FindMarkers(seurat_object,
                    ident.1 = 'KD',
                    ident.2 = 'Wt',
                    group.by = 'genotype',
                    subset.ident = cell_type) |>
    rownames_to_column('gene') |>
    left_join(gene_name_to_entrez_id, by = 'gene')

  output_subfolder_kd_cell_type <- file.path(output_folder, 'KD_effect_by_cell_type', cell_type)
  dir.create(output_subfolder_kd_cell_type, recursive = TRUE)

  gene_sct_count_total <- rowSums(subset(seurat_object, idents = cell_type)@assays$
                                    SCT$
                                    counts)
  expressed_gene_names_in_cell_type <- names(gene_sct_count_total[gene_sct_count_total > 0])


  save_de_and_go_enrichment(de, expressed_gene_names_in_cell_type, output_subfolder_kd_cell_type)

  # Cell type vs all other cells
  if ((!'Wt' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type])) |
    (!'KD' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type]))) {
    print(paste(cell_type, 'is not present in both genotype groups!'))
    next
  }
  de_cell_type_vs_others <- FindMarkers(seurat_object,
                                        ident.1 = cell_type) |>
    rownames_to_column('gene') |>
    left_join(gene_name_to_entrez_id, by = 'gene')

  output_subfolder_cell_type_vs_others <- file.path(output_folder, 'cell_type_vs_others', cell_type)
  dir.create(output_subfolder_cell_type_vs_others, recursive = TRUE)

  save_de_and_go_enrichment(de_cell_type_vs_others,
                            expressed_gene_names_all_cells,
                            output_subfolder_cell_type_vs_others)

}


