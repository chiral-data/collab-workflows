FROM condaforge/miniforge3:latest

ENV USER=silva HOME=/tmp MPLCONFIGDIR=/tmp/matplotlib

# AmberTools first — pins ICU and HDF5 versions
RUN conda install -y -c conda-forge \
        python=3.11 ambertools && conda clean -afy

# RDKit + packmol after AmberTools to respect ICU pin
RUN conda install -y -c conda-forge \
        rdkit packmol && conda clean -afy

# acpype for Amber → GROMACS topology conversion
RUN pip install --no-cache-dir acpype

WORKDIR /workspace
