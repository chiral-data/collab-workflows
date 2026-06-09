FROM continuumio/miniconda3:latest

# Runtime dependencies needed by the GNINA prebuilt binary
RUN apt-get update && apt-get install -y \
    wget \
    libboost-filesystem1.74.0 \
    libboost-program-options1.74.0 \
    libboost-system1.74.0 \
    libboost-thread1.74.0 \
    libboost-serialization1.74.0 \
    libboost-iostreams1.74.0 \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN conda install -y -c conda-forge \
    python=3.11 \
    rdkit \
    numpy \
    pandas \
    && conda clean -afy

# GNINA prebuilt binary (v1.3, Linux x86_64)
# Falls back to CPU automatically when no GPU is present
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina \
    -O /usr/local/bin/gnina && chmod +x /usr/local/bin/gnina

WORKDIR /workspace
