# Use NVIDIA Pre-built Gromacs Image 
FROM nvcr.io/hpc/gromacs:2023.2

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends zip \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a working directory
WORKDIR /workspace

ENV PATH="/usr/local/gromacs/avx2_256/bin:${PATH}"

# Default command, used for verification
CMD ["gmx", "--version"]
