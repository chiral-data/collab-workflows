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

# Switch to root to perform file operations
USER root

# Copy source code (respecting .dockerignore if present, but we copy widely here)
# We copy specific directories to avoid clutter
COPY 1_data_prep /app/1_data_prep
COPY 2_model_train /app/2_model_train
COPY 3_analyze_overfitting /app/3_analyze_overfitting
COPY 4_predict_from_csv /app/4_predict_from_csv
COPY input_files /app/input_files

# Create output directories to ensure they exist/persist
RUN mkdir -p /app/1_data_prep/outputs \
    /app/2_model_train/outputs \
    /app/3_analyze_overfitting/outputs \
    /app/4_predict_from_csv/outputs

# Make run scripts executable
RUN chmod +x /app/1_data_prep/run.sh \
    /app/2_model_train/run.sh \
    /app/3_analyze_overfitting/run.sh \
    /app/4_predict_from_csv/run.sh

# NOTE: Keep running as root so mounted volumes from host are writable.
# Proper fix: silva executor should set container user to match host UID:GID.

# Expose data volumes for explicit input/output handling
VOLUME ["/app/input_files", "/app/4_predict_from_csv/inputs"]

# Default command: Open a shell
CMD ["/bin/bash"]
