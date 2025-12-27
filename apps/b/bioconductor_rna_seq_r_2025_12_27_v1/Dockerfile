# RNA-seq analysis Docker image
# Based on Bioconductor with DESeq2, Rsubread, and visualization packages
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages
RUN Rscript -e '\
    BiocManager::install(c( \
        "Rsubread", \
        "DESeq2", \
        "ashr" \
    ), ask = FALSE, update = FALSE)'

# Install CRAN packages
RUN Rscript -e '\
    install.packages(c( \
        "tidyverse", \
        "ggrepel", \
        "ggrastr" \
    ), repos = "https://cloud.r-project.org")'

# Set working directory
WORKDIR /workspace

# Verify installation
RUN Rscript -e 'library(Rsubread); library(DESeq2); library(ashr); library(tidyverse); library(ggrepel); library(ggrastr); cat("All packages loaded successfully\n")'

# Default command
CMD ["R", "--version"]
