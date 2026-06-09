FROM continuumio/miniconda3:latest

RUN conda install -y -c conda-forge \
    python=3.11 \
    pandas \
    scipy \
    numpy \
    openbabel \
    && conda clean -afy

WORKDIR /workspace
