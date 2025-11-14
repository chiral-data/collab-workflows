# Use Ubuntu 22.04 base image to match project pattern
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies, python3, and jq for JSON config parsing
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    python3 \
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba (lightweight conda alternative)
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
SHELL ["/bin/bash", "-l", "-c"]

RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin/ --strip-components=1 bin/micromamba

# Initialize micromamba and install RNAhybrid in one step
RUN micromamba create -y -n base -c defaults -c bioconda -c conda-forge rnahybrid=2.1.2 && \
    micromamba clean --all --yes

# Set up environment activation
ENV MAMBA_ACTIVATE_PREFIX=$MAMBA_ROOT_PREFIX
RUN echo "eval \"\$(micromamba shell hook --shell bash)\"" >> ~/.bashrc && \
    echo "micromamba activate base" >> ~/.bashrc

# Create a working directory
WORKDIR /workspace

# Default command to verify installation
CMD ["RNAhybrid", "-h"]
