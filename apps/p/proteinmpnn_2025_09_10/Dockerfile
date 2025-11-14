# Use NVIDIA CUDA base image with Python
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a symbolic link for python
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Install PyTorch with CUDA support and NumPy
# Installing torch only first (without torchvision/torchaudio to reduce size)
RUN pip install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cu126
RUN pip install --no-cache-dir numpy

# Clone ProteinMPNN repository
WORKDIR /opt
RUN git clone https://github.com/dauparas/ProteinMPNN.git

# Add ProteinMPNN to Python path
ENV PYTHONPATH="/opt/ProteinMPNN:${PYTHONPATH}"

# Create a working directory
WORKDIR /workspace

# Default command, used for verification
CMD ["python", "/opt/ProteinMPNN/protein_mpnn_run.py", "--help"]
