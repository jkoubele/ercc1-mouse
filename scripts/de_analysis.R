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

gene_sct_count_total <- rowSums(seurat_object@assays$SCT$counts)
expressed_gene_names <- names(gene_sct_count_total[gene_sct_count_total > 0])
expressed_gene_ids <-  bitr(expressed_gene_names,
                                fromType = "SYMBOL",
                                toType = "ENTREZID", 
                                OrgDb = org.Mm.eg.db)


de_genotype <- FindMarkers(seurat_object,
                           ident.1 = 'KD',
                           ident.2 = 'Wt',
                           group.by = 'genotype') |>
  rownames_to_column('gene')

de_genotype_up <- de_genotype |> dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05)
de_genotype_down <- de_genotype |> dplyr::filter(avg_log2FC < 0, p_val_adj < 0.05)

de_genotype_up_gene_id <-  bitr(de_genotype_up$gene,
                                fromType = "SYMBOL",
                                toType = "ENTREZID", 
                                OrgDb = org.Mm.eg.db)

de_genotype_down_gene_id <-  bitr(de_genotype_down$gene,
                                fromType = "SYMBOL",
                                toType = "ENTREZID", 
                                OrgDb = org.Mm.eg.db)

go_mf_up <- enrichGO(gene = de_genotype_up_gene_id$ENTREZID, 
                         OrgDb = org.Mm.eg.db, 
                         ont = "MF",
                        universe = expressed_gene_ids$ENTREZID,
                         readable = TRUE)

go_mf_down <- enrichGO(gene = de_genotype_down_gene_id$ENTREZID, 
                     OrgDb = org.Mm.eg.db, 
                     ont = "MF",
                     universe = expressed_gene_ids$ENTREZID,
                     readable = TRUE)

barplot_go_mf_up <- barplot(go_mf_up, showCategory=10, title='GO enrichment (MF)')
barplot_go_mf_down <- barplot(go_mf_down, showCategory=10, title='GO enrichment (MF)')

dotplot_go_mf_up <- barplot(go_mf_up, showCategory=10, title='GO enrichment (MF)')



write_tsv(de_genotype, file.path(output_folder, "de_genotype_all.tsv"))
write_tsv(de_genotype_up, file.path(output_folder, "de_genotype_up.tsv"))
write_tsv(de_genotype_down, file.path(output_folder, "de_genotype_down.tsv"))

Idents(seurat_object) <- "cell_type"
for (cell_type in unique(seurat_object$cell_type)) {
  if((! 'Wt' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type])) |
     (! 'KD' %in% unique(seurat_object$genotype[seurat_object$cell_type == cell_type])) ){
    next # cell type is not present in both genotype groups
  }

  de_same_cell_type <- FindMarkers(seurat_object,
                                   ident.1 = 'KD',
                                   ident.2 = 'Wt',
                                   group.by = 'genotype',
                                   subset.ident = cell_type) |> 
    rownames_to_column('gene')

  de_same_cell_type_up <- de_same_cell_type |> dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05)
  de_same_cell_type_down <- de_same_cell_type |> dplyr::filter(avg_log2FC < 0, p_val_adj < 0.05)
  
  output_subfolder <- file.path(output_folder, 'de_within_cell_type', cell_type)
  dir.create(output_subfolder, recursive=TRUE)
  write_tsv(de_same_cell_type, file.path(output_subfolder, "de_all.tsv"))
  write_tsv(de_same_cell_type_up, file.path(output_subfolder, "de_up.tsv"))
  write_tsv(de_same_cell_type_down, file.path(output_subfolder, "de_down.tsv"))
  
  # TO-DO: replace with FindAllMarkers() ?
  de_cell_type_vs_others <- FindMarkers(seurat_object,
                                   ident.1 = cell_type) |>
    rownames_to_column('gene')
  de_cell_type_vs_others_up <- de_cell_type_vs_others |> dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05)
  de_cell_type_vs_others_down <- de_cell_type_vs_others |> dplyr::filter(avg_log2FC < 0, p_val_adj < 0.05)

  output_subfolder_cell_type_vs_others <- file.path(output_folder, 'de_cell_type_vs_others', cell_type)
  dir.create(output_subfolder_cell_type_vs_others, recursive=TRUE)

  write_tsv(de_cell_type_vs_others, file.path(output_subfolder_cell_type_vs_others, "de_all.tsv"))
  write_tsv(de_cell_type_vs_others_up, file.path(output_subfolder_cell_type_vs_others, "de_up.tsv"))
  write_tsv(de_cell_type_vs_others_down, file.path(output_subfolder_cell_type_vs_others, "de_down.tsv"))
  
}


# Idents(seurat_object) <- "intervention_cell_type"
# 
# de_1 <- FindMarkers(seurat_object, ident.1 = "Hematopoietic_cells", ident.2 = "Proximal_tubule_cells")
# 
# cell_mapping <- data.frame(
#   simplified = colnames(seurat_object),
#   barcode = colnames(GetAssayData(seurat_object, assay = "SCT", layer = "data")),
#   stringsAsFactors = FALSE
# )group_1_barcodes <- cell_mapping$barcode[cell_mapping$simplified %in% WhichCells(seurat_object, idents = "Wt_Distal_tubule_cells")]
# group_2_barcodes <- cell_mapping$barcode[cell_mapping$simplified %in% WhichCells(seurat_object, idents = "KD_Distal_tubule_cells")]sct_assay <- seurat_object@assays$SCTde_1 <- FindMarkers(
#   object = sct_assay,
#   cells.1 = group_1_barcodes,
#   cells.2 = group_2_barcodes,
#   test.use = "wilcox"
# )de_1_up <- de_1 |> dplyr::filter(avg_log2FC > 0, p_val_adj<0.05)
# de_1_down <- de_1 |> dplyr::filter(avg_log2FC < 0, p_val_adj<0.05)
