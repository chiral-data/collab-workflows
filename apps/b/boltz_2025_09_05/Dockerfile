# Use NVIDIA CUDA base image with Python
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    BOLTZ_CACHE="/opt/boltz_cache"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    git \
    wget \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a symbolic link for python
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Create boltz cache directory
RUN mkdir -p ${BOLTZ_CACHE}

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Method 1: Install Boltz first, then fix dependencies
RUN pip install --no-cache-dir boltz

# Force correct versions AFTER Boltz installation
RUN pip install --no-cache-dir --force-reinstall \
    "numpy>=2.0,<2.2" \
    "numba>=0.59.0"

# Method 2 (Alternative): Use constraints file
# Create constraints file
# RUN echo "numpy<2.2" > /tmp/constraints.txt
# RUN pip install --no-cache-dir -c /tmp/constraints.txt boltz

# Method 3 (Alternative): Install without deps then manually install
# RUN pip install --no-cache-dir --no-deps boltz
# RUN pip install --no-cache-dir \
#     "numpy>=2.0,<2.2" \
#     "numba>=0.59.0" \
#     torch torchvision \
#     pytorch-lightning \
#     # ... other dependencies

# Install additional packages
RUN pip install --no-cache-dir tensorboard pytorch-lightning==2.5.2

# Verify NumPy version
RUN python -c "import numpy; print(f'NumPy version: {numpy.__version__}')"

# Pre-download Boltz-2 model weights to cache directory
RUN cd ${BOLTZ_CACHE} && \
    wget -q https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt && \
    wget -q https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt && \
    wget -q https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar && \
    tar -xf mols.tar

# Create a working directory
WORKDIR /workspace

# Default command, used for verification
CMD ["boltz", "--help"]
