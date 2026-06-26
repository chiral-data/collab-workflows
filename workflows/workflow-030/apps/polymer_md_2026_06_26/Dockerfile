FROM continuumio/miniconda3

RUN conda install -y python=3.11 setuptools && conda clean -afy

ENV USER=silva \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

# RDKit (oligomer build) + AmberTools (GAFF2 charges) + packmol (cell packing) + acpype (GROMACS conversion)
RUN conda install -y -c conda-forge \
        rdkit \
        ambertools \
        packmol \
        acpype \
    && conda clean -afy

WORKDIR /workspace
