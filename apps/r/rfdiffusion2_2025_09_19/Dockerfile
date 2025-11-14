# Use NVIDIA CUDA base image for GPU support
FROM nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      git \
      python3 \
      python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Clone RFdiffusion2 repository to /opt/RFdiffusion2
WORKDIR /opt
RUN git clone https://github.com/RosettaCommons/RFdiffusion2.git RFdiffusion2

# Set working directory
WORKDIR /opt/RFdiffusion2

# Download model weights
RUN python3 setup.py

# Set environment
ENV PYTHONPATH=/opt/RFdiffusion2

# Default command shows help
CMD ["python3", "rf_diffusion/benchmark/pipeline.py", "--help"]