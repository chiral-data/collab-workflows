# Dockerfile to build a custom image for ViennaRNA 
#
# 2025-11-15: add ghostscript, pdf2svg to convert .ps files (Qin Wan rogerwq@gmai.com)
#

from ubuntu:22.04

# set environment variables
env debian_frontend=noninteractive \
    pythondontwritebytecode=1 \
    pythonunbuffered=1

# install system dependencies and wget for micromamba
run apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    jq \
    ghostscript \
    pdf2svg \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba (lightweight conda alternative)
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
SHELL ["/bin/bash", "-l", "-c"]

RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin/ --strip-components=1 bin/micromamba

# Initialize micromamba and install ViennaRNA in one step
RUN micromamba create -y -n base -c defaults -c bioconda -c conda-forge viennarna=2.7.0 && \
    micromamba clean --all --yes

# Set up environment activation
ENV MAMBA_ACTIVATE_PREFIX=$MAMBA_ROOT_PREFIX
RUN echo "eval \"\$(micromamba shell hook --shell bash)\"" >> ~/.bashrc && \
    echo "micromamba activate base" >> ~/.bashrc

# Create a working directory
WORKDIR /workspace

# Default command to verify installation
CMD ["RNAfold", "--version"]
