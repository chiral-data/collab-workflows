FROM condaforge/miniforge3:latest

ENV USER=silva \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# AmberTools, RDKit, packmol in one solve to let libmamba resolve conflicts —
# separate RUN layers re-solve the whole env each time and can silently
# unlink an earlier package with no error (see #212)
RUN mamba install -y -c conda-forge \
        python=3.11 \
        ambertools \
        rdkit \
        packmol \
    && mamba clean -afy

# acpype via pip (thin pure-Python wrapper around AmberTools)
RUN pip install --no-cache-dir acpype

WORKDIR /workspace
