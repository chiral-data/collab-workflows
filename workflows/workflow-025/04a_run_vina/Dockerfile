FROM continuumio/miniconda3:latest

RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

RUN mamba install -y -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    rdkit \
    openbabel \
    vina \
    && mamba clean -afy

WORKDIR /workspace
