FROM continuumio/miniconda3:latest

RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

RUN mamba install -y -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    pdbfixer \
    openmm \
    rdkit \
    openbabel \
    xtb \
    && mamba clean -afy

RUN pip install --no-cache-dir biopython meeko

WORKDIR /workspace
