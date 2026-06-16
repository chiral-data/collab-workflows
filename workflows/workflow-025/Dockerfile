FROM condaforge/miniforge3:latest

# CLI tools and Boost runtime libs for GNINA binary
RUN apt-get update && apt-get install -y \
    wget \
    xz-utils \
    openbabel \
    libboost-filesystem-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-iostreams-dev \
    libnvtoolsext1 \
    && rm -rf /var/lib/apt/lists/*

# xtb and openbabel are CLI tools installed outside conda to avoid solver conflicts
RUN mamba install -y -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    scipy \
    pdbfixer \
    openmm \
    rdkit \
    && mamba clean -afy

RUN pip install --no-cache-dir biopython meeko \
    nvidia-cuda-runtime-cu12 \
    nvidia-cublas-cu12 \
    nvidia-cusparse-cu12 \
    nvidia-cufft-cu12 \
    nvidia-cusolver-cu12 \
    nvidia-nvtx-cu12

# GNINA dynamically links CUDA libraries even in CPU mode; expose them all
ENV LD_LIBRARY_PATH=/opt/conda/lib/python3.11/site-packages/nvidia/cuda_runtime/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cublas/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cusparse/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cufft/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cusolver/lib:/opt/conda/lib/python3.11/site-packages/nvidia/nvtx/lib

# xtb v6.7.1 prebuilt binary
RUN wget -q https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz \
    -O /tmp/xtb.tar.xz \
    && tar -xf /tmp/xtb.tar.xz -C /tmp \
    && mv /tmp/xtb-*/bin/xtb /usr/local/bin/xtb \
    && rm -rf /tmp/xtb*

# AutoDock Vina v1.2.5 prebuilt binary
RUN wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    -O /usr/local/bin/vina && chmod +x /usr/local/bin/vina

# GNINA v1.3 prebuilt binary (falls back to CPU without GPU)
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina \
    -O /usr/local/bin/gnina && chmod +x /usr/local/bin/gnina

WORKDIR /workspace
