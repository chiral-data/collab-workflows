FROM continuumio/miniconda3:26.7.1-1

# Runtime dependencies needed by the GNINA prebuilt binary
RUN apt-get update && apt-get install -y \
    wget \
    libboost-filesystem-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-iostreams-dev \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN conda install -y -c conda-forge \
    python=3.11 \
    rdkit \
    numpy \
    pandas \
    && conda clean -afy

RUN pip install --no-cache-dir \
    nvidia-cuda-runtime-cu12 \
    nvidia-cublas-cu12 \
    nvidia-cusparse-cu12 \
    nvidia-cufft-cu12 \
    nvidia-cusolver-cu12 \
    nvidia-nvtx-cu12==12.8.90  # 12.9+ dropped libnvToolsExt.so.1, which gnina links

# GNINA dynamically links CUDA libraries even in CPU mode; expose them all
ENV LD_LIBRARY_PATH=/opt/conda/lib/python3.11/site-packages/nvidia/cuda_runtime/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cublas/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cusparse/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cufft/lib:/opt/conda/lib/python3.11/site-packages/nvidia/cusolver/lib:/opt/conda/lib/python3.11/site-packages/nvidia/nvtx/lib

# GNINA prebuilt binary (v1.3, Linux x86_64)
# Falls back to CPU automatically when no GPU is present
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina \
    -O /usr/local/bin/gnina && chmod +x /usr/local/bin/gnina

WORKDIR /workspace
