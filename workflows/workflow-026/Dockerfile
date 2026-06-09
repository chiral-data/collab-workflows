FROM ghcr.io/sokrypton/colabfold:1.6.1-cuda12

# PyMOL for structure visualization (Node 04)
RUN mamba install -y -c conda-forge pymol-open-source=3.1.0 && \
    mamba clean -afy

# Download AlphaFold2 weights into /cache so they're baked into the image
RUN python -m colabfold.download
