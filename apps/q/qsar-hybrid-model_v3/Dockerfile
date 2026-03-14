FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /app

# Metadata
LABEL version="3.0"
LABEL description="QSAR Hybrid Model V3 - Interactive Dashboards & Prediction Pipeline"
LABEL maintainer="QSAR Abdo"

# Copy environment file
COPY final_environment.yml /tmp/env.yaml

# Install dependencies
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Activate environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cudnn/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cublas/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cuda_runtime/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cuda_cupti/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cuda_nvcc/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cufft/lib:/opt/conda/lib/python3.10/site-packages/nvidia/curand/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cusolver/lib:/opt/conda/lib/python3.10/site-packages/nvidia/cusparse/lib:/opt/conda/lib/python3.10/site-packages/nvidia/nccl/lib:/opt/conda/lib/python3.10/site-packages/nvidia/nvjitlink/lib"

# TensorFlow & NVIDIA Optimizations
ENV TF_ENABLE_ONEDNN_OPTS=0
ENV TF_CPP_MIN_LOG_LEVEL=2
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

WORKDIR /workspace

CMD ["/bin/bash"]
