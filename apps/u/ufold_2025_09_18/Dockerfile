# Use NVIDIA CUDA base image compatible with UFold requirements
FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies including OpenCV requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3-pip \
    git \
    wget \
    build-essential \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create symbolic links for python
RUN ln -sf /usr/bin/python3.11 /usr/bin/python3
RUN ln -sf /usr/bin/python3.11 /usr/bin/python

# Upgrade pip
RUN python -m pip install --no-cache-dir --upgrade pip

# Install PyTorch 2.0.1 with CUDA 11.8 support (matching UFold requirements)
RUN pip install --no-cache-dir \
    torch==2.0.1 \
    torchvision==0.15.2 \
    --index-url https://download.pytorch.org/whl/cu118

# Install UFold dependencies with exact versions from UFold.yaml
RUN pip install --no-cache-dir \
    numpy==1.24.3 \
    pandas==2.1.1 \
    scipy==1.11.3 \
    opencv-python-headless==4.8.1.78 \
    matplotlib==3.8.0 \
    scikit-learn==1.3.0 \
    huggingface-hub==0.17.3 \
    timm==0.9.7 \
    munch

# Clone UFold repository
WORKDIR /opt
RUN git clone https://github.com/uci-cbcl/UFold.git

# Copy data and models folders into the UFold directory
COPY data/ /opt/UFold/data/
COPY models/ /opt/UFold/models/

# Add UFold to Python path
ENV PYTHONPATH="/opt/UFold:${PYTHONPATH}"

# Create a working directory
WORKDIR /workspace

# Default command, used for verification
CMD ["python", "/opt/UFold/ufold_predict.py", "--help"]