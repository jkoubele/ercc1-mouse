# ercc1-mouse

Reads were aligned to the mm10 genome (aka GRCm38) using Cellranger, via the script [batch_align_cellranger.sh](cluster_pipeline/batch_align_cellranger.sh).
I used GRCm38 instead of the newer GRCm39, since I intend to use [SComatic](https://github.com/cortes-ciriano-lab/SComatic) in the downstream analysis, and SComatic uses RNA editing database which seems to be available only for GRCm38.

The standard single-cell data analysis is implemented in R, mostly using the Seurat package. Running the R scripts can 
be done from the customized RStudio docker image, which can be built by:
```
docker build -t custom_bioconductor ./docker_files/bioconductor/
```

The RStudio server can be then run by:
```
 docker run -it -e PASSWORD=pass -e USERID=$(id -u) -e GROUPID=$(id -g) -p 8787:8787 -v /cellfile/datapublic/jkoubele/ercc1:/ercc1 custom_bioconductor
```

The analysis pipeline consists of several steps, each run by separate script. The results of these
steps are saved to disk (serializing the Seurat objects to .rds files).

1. **Ambient RNA removal** using [decontX](https://www.bioconductor.org/packages/release/bioc/html/decontX.html) can be run using the script [ambient_rna_removal.R](scripts/ambient_rna_removal.R).
2. **Doublets detection** using [scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) can be run by using the script [doublet_detection.R](scripts/doublet_detection.R).
3. **Cell QC**: filter out outlier cells based on the number of detected genes, UMIs and mitochondrial content, using the script [cell_quality_control.R](scripts/cell_quality_control.R).
4. **Clustering and dim. reduction** is performed by the script [clustering.R](scripts/clustering.R).