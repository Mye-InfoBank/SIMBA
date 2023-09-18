FROM rocker/r2u:latest

RUN Rscript -e 'install.packages(c("BiocManager", "Seurat", "anndata", "MatrixExtra", "devtools"))'

RUN Rscript -e 'BiocManager::install(c("celda", "SingleR"))'

RUN Rscript -e 'devtools::install_github("https://github.com/mimi3421/ROGUE/tree/patch-1")'

RUN apt install -y git libffi-dev libssl-dev libsqlite3-dev

RUN Rscript -e 'reticulate::install_python()'

RUN Rscript -e 'reticulate::virtualenv_create()'

RUN Rscript -e 'anndata::install_anndata()'