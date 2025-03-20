FROM bioconductor/bioconductor_docker
RUN R -e "install.packages(c('tidyverse', 'rjson', 'argparse', 'HGNChelper', 'openxlsx', 'homologene', 'Polychrome'))"
RUN R -e "BiocManager::install(c('decontX', 'DropletUtils', 'scDblFinder', 'glmGamPoi', 'biomaRt'),  \
    dependencies = TRUE, ask = FALSE)"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"