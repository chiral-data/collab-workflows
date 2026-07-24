FROM condaforge/miniforge3:latest

ENV USER=silva \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

# AmberTools first — it pins icu and other system libs that rdkit must respect
RUN mamba install -y -c conda-forge \
        python=3.11 \
        ambertools \
    && mamba clean -afy

# RDKit + packmol into the same environment, solved after AmberTools pins
RUN mamba install -y -c conda-forge \
        rdkit \
        packmol \
    && mamba clean -afy

# acpype via pip (thin pure-Python wrapper around AmberTools)
RUN pip install --no-cache-dir acpype

WORKDIR /workspace
