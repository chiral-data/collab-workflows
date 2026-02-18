# BoltzGen - AI-based protein binder design
# Based on tested Dockerfile from container-images-for-potter
# Use NVIDIA CUDA devel image (needed for cuequivariance compilation)
FROM nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_EXTRA_INDEX_URL=https://download.pytorch.org/whl/cu121 \
    HF_HOME="/opt/boltzgen_cache"

# Install system dependencies and Python 3.11
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3.11-distutils \
    python3-pip \
    build-essential \
    git \
    jq \
    libgl1 \
    libhdf5-dev \
    libboost-all-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create symbolic links for python3.11
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# Upgrade pip
RUN python3.11 -m pip install --upgrade pip setuptools wheel

# Create cache directory for model weights
RUN mkdir -p ${HF_HOME}

# Install BoltzGen and dependencies (PyTorch resolved via cu121 extra index)
RUN pip install --no-cache-dir boltzgen pyyaml

# Verify installation
RUN python -c "import boltzgen; print('BoltzGen installed successfully')"

# Pre-download model weights and datasets to cache directory
RUN boltzgen download --cache ${HF_HOME} all

# Create working directory
WORKDIR /workspace

# Keep container running for silva exec
CMD ["tail", "-f", "/dev/null"]
