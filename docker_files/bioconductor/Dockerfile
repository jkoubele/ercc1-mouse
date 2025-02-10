FROM bioconductor/bioconductor_docker
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('rjson')"
RUN Rscript -e 'install.packages("argparse")'
RUN Rscript -e "BiocManager::install('decontX', dependencies = TRUE, ask = FALSE)"
RUN Rscript -e "BiocManager::install('DropletUtils', dependencies = TRUE, ask = FALSE)"
RUN Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN Rscript -e 'BiocManager::install("scDblFinder", dependencies = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("glmGamPoi", dependencies = TRUE, ask = FALSE)'
RUN Rscript -e 'install.packages("HGNChelper")'
RUN Rscript -e 'install.packages("openxlsx")'
