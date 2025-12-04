FROM nvidia/cuda:12.6.0-cudnn-devel-ubuntu24.04

WORKDIR /root

# Install build tools, libraries, and Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libboost-all-dev \
    libeigen3-dev \
    libopenbabel-dev \
    python3-dev \
    python3-setuptools \
    wget \
    unzip \
    python3-numpy \
    python3-scipy \
    cython3 \
    python3-openbabel \
    libprotobuf-dev \
    protobuf-compiler \
    libjsoncpp-dev \
    python3-pytest

# Install a newer version of CMake from Kitware APT repository
RUN apt-get remove -y --purge --auto-remove cmake && \
    apt-get update && \
    apt-get install -y ca-certificates gpg wget && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null && \
    echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ noble main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null && \
    apt-get update && \
    apt-get install -y cmake

# Clone and build Open Babel (if not already provided by libopenbabel-dev)
RUN git clone --depth 1 --branch v1.3 https://github.com/gnina/gnina.git
WORKDIR /root/gnina
RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

# Set environment variables for gnina
ENV PATH="/root/gnina/build/bin:${PATH}"

# Create a working directory
WORKDIR /workspace

# Command to show gnina version
CMD ["gnina", "--version"]
