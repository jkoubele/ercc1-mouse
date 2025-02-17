library(tidyverse)


path_prefix <- "/ercc1/data"
output_folder <- file.path(path_prefix, "CCG_data")


sample_annotation <- data.frame(
  CCG_sample_id = paste0("A006850389_23250", 3:6),
  original_sample_name = paste0("24T256-", 1:4),
  sample_name = c("KD_rep_1", "KD_rep_2", "Wt_rep_1", "Wt_rep_2"),
  genotype = c("KD", "KD", "Wt", "Wt")
)

# Save to TSV file
write_tsv(sample_annotation, file.path(output_folder, "sample_annotation.tsv"))