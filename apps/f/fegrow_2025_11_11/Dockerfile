# Ubuntu-based Dockerfile for workflow-functional-group environment
FROM ubuntu:24.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    libgl1-mesa-dri \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    libgomp1 \
    libgcc-s1 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Initialize conda
RUN conda init bash

# Accept conda terms of service and configure channels
RUN conda config --set channel_priority strict && \
    conda config --add channels conda-forge && \
    conda config --set auto_activate_base false && \
    conda config --set always_yes true

# Accept terms of service for required channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install mamba for faster package management
RUN conda install -c conda-forge mamba -y

# Clone the project
RUN git clone https://github.com/cole-group/FEgrow.git /opt/FEgrow

# Create conda environment from environment.yml
RUN mamba env create -f /opt/FEgrow/environment.yml

# Install essential packages first
RUN mamba run -n fegrow pip install numpy pandas matplotlib seaborn scipy scikit-learn biopython

# Try to install additional packages from requirements.txt (skip problematic ones)
RUN mamba run -n fegrow pip install -r /opt/FEgrow/requirements.txt --ignore-requires-python || echo "Some packages failed to install, continuing..."

# Install FEgrow package 
RUN mamba run -n fegrow pip install --no-deps /opt/FEgrow

# Create a wrapper script to activate the environment
RUN echo '#!/bin/bash\nsource /opt/conda/etc/profile.d/conda.sh\nconda activate fegrow\nexec "$@"' > /usr/local/bin/workflow-run && \
    chmod +x /usr/local/bin/workflow-run

# Set default conda environment
RUN echo 'conda activate fegrow' >> /root/.bashrc

# Set the default command
CMD ["/bin/bash"]

# Expose port for Jupyter notebook if needed
EXPOSE 8888

# Add labels for metadata
LABEL maintainer="Workflow Functional Group Docker Environment"
LABEL description="Ubuntu-based Docker environment for workflow-functional-group"

# Set working directory
WORKDIR /workspace
