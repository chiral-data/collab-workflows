FROM condaforge/miniforge3:latest

ENV USER=silva HOME=/tmp MPLCONFIGDIR=/tmp/matplotlib

# AmberTools first — pins ICU and HDF5 versions
RUN micromamba install -y -n base -c conda-forge \
        python=3.11 ambertools && micromamba clean -afy

# RDKit + packmol after AmberTools to respect ICU pin
RUN micromamba install -y -n base -c conda-forge \
        rdkit packmol && micromamba clean -afy

# acpype for Amber → GROMACS topology conversion
RUN pip install --no-cache-dir acpype

WORKDIR /workspace
