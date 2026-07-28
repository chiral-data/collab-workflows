FROM mambaorg/micromamba:latest
USER root

ENV HOME=/tmp MPLCONFIGDIR=/tmp/matplotlib \
    PATH=/opt/conda/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# AmberTools, RDKit, packmol in one solve to let libmamba resolve conflicts
RUN micromamba install -y -n base -c conda-forge \
        python=3.11 ambertools rdkit packmol && \
    micromamba clean -afy

# acpype for Amber → GROMACS topology conversion
RUN /opt/conda/bin/pip install --no-cache-dir acpype

WORKDIR /workspace
