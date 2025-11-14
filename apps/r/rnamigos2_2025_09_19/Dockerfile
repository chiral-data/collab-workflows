# Use Python 3.11 base image (Ubuntu 22.04 based for compatibility)
FROM python:3.11-slim-bookworm

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    build-essential \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    libgomp1 \
    libxrender1 \
    libxinerama1 \
    libxi6 \
    libfontconfig1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Install numpy first (required version)
RUN pip install --no-cache-dir numpy==1.26

# Install PyTorch CPU version and related packages with specific versions
RUN pip install --no-cache-dir \
    torch==2.2.2+cpu \
    torchaudio==2.2.2 \
    torchdata==0.7.1 \
    torchvision==0.17.2 \
    --index-url https://download.pytorch.org/whl/cpu

# Install DGL from specific wheel repository
RUN pip install --no-cache-dir dgl -f https://data.dgl.ai/wheels/torch-2.2/repo.html

# Clone RNAmigos2 repository
WORKDIR /opt
RUN git clone https://github.com/cgoliver/rnamigos2.git

# Install RNAmigos2 requirements
WORKDIR /opt/rnamigos2
RUN pip install --no-cache-dir -r requirements.txt

# Add RNAmigos2 to Python path
ENV PYTHONPATH="/opt/rnamigos2:$PYTHONPATH"

# Create a working directory
WORKDIR /workspace

# Copy sample files for testing (if needed)
RUN cp -r /opt/rnamigos2/data/sample_files /workspace/sample_files

# Pre-download the RNA-FM pretrained model to avoid downloading it on every run
RUN mkdir -p /root/.cache/torch/hub/checkpoints && \
    wget -q --show-progress --progress=bar:force:noscroll \
    "https://proj.cse.cuhk.edu.hk/rnafm/api/download?filename=RNA-FM_pretrained.pth" \
    -O /root/.cache/torch/hub/checkpoints/RNA-FM_pretrained.pth && \
    echo "RNA-FM model downloaded successfully"

# Default command to show help
CMD ["python", "/opt/rnamigos2/rnamigos/inference.py", "--help"]