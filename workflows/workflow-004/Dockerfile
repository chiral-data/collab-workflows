# Unified Dockerfile for AutoDock Vina Workflow
# This single image contains all tools needed for all 8 nodes

FROM continuumio/miniconda3:latest

# Update conda and install mamba for faster package resolution
RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

# Install Python and core scientific packages
RUN mamba install -y -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    openpyxl \
    && mamba clean -afy

# Install molecular modeling packages
RUN mamba install -y -c conda-forge \
    pdbfixer \
    openmm \
    rdkit \
    openbabel \
    && mamba clean -afy

# Install AutoDock Vina
RUN mamba install -y -c conda-forge vina && \
    mamba clean -afy

# Install Python packages via pip
RUN pip install --no-cache-dir \
    biopython \
    meeko \
    gemmi

# Install Java for P2Rank
RUN apt-get update && \
    apt-get install -y default-jre-headless wget && \
    rm -rf /var/lib/apt/lists/*

# Download and install P2Rank
RUN wget https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz && \
    tar -xzf p2rank_2.4.2.tar.gz && \
    mv p2rank_2.4.2 /opt/p2rank && \
    rm p2rank_2.4.2.tar.gz

ENV PATH="/opt/p2rank:${PATH}"

WORKDIR /workspace

# Metadata
LABEL description="Unified image for AutoDock Vina workflow - all 8 nodes"
LABEL version="1.0"
LABEL tools="Biopython, PDBFixer, OpenMM, Meeko, RDKit, OpenBabel, Vina, P2Rank"
