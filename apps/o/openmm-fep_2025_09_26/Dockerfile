# Use NVIDIA CUDA base image with Python
FROM nvidia/cuda:12.2.2-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (Conda) at /miniforge3
ENV CONDA_DIR=/miniforge3
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm -f /tmp/miniforge.sh

# Put conda on PATH
ENV PATH=${CONDA_DIR}/bin:${PATH}

# Improve conda robustness and cleanup
RUN conda config --set channel_priority strict && \
    conda config --set always_yes yes && \
    conda update -n base -c conda-forge conda && \
    conda clean -afy

# Install OpenMM FEP toolkits and dependencies
RUN conda install -c conda-forge \
    openmm \
    cuda-version=12.2 \
    openmmtools \
    openff-toolkit \
    openmmforcefields \
    openbabel \
    ambertools \
    pdbfixer \
    pymbar \
    pandas \
    tqdm \
    hydra-core \
    omegaconf

# Set AMBERHOME environment variable for AmberTools
ENV AMBERHOME=/miniforge3

# Source amber.sh in bashrc for proper AmberTools setup
RUN echo "source ${AMBERHOME}/amber.sh" >> ~/.bashrc

# Verify OpenMM installation
RUN python -m openmm.testInstallation

# Create a working directory
WORKDIR /workspace

# Default command, used for verification
CMD ["python", "-m", "openmm.testInstallation"]