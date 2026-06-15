# ============================================================
# workflow-010-antibody-structure-prediction
# CUDA 12.1 image with ProtT5 support for ABB3-LM
# ============================================================

FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# ---- System dependencies -----------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.10 \
        python3-pip \
        build-essential \
        git \
        wget \
        curl \
        dos2unix \
    && ln -sf python3.10 /usr/bin/python \
    && ln -sf pip3 /usr/bin/pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workflow

# ---- Python dependencies -----------------------------------------------
RUN pip install --no-cache-dir --upgrade pip

# PyTorch with CUDA 12.1 support
RUN pip install --no-cache-dir \
    torch==2.4.0 \
    --index-url https://download.pytorch.org/whl/cu121

# Core ML / bio dependencies
RUN pip install --no-cache-dir \
    numpy \
    scipy \
    biopython \
    ml-collections \
    einops

# ProtT5 dependencies (transformers + sentencepiece for Node 02 / ABB3-LM)
RUN pip install --no-cache-dir \
    transformers \
    sentencepiece

# pdbfixer + openmm (required by abodybuilder3's openfold relax module)
RUN pip install --no-cache-dir \
    pdbfixer \
    openmm

# ---- ABodyBuilder3 -------------------------------------------------------
RUN git clone https://github.com/Exscientia/abodybuilder3.git abodybuilder3 \
    && cd abodybuilder3 \
    && mkdir -p output \
    && wget -q https://zenodo.org/records/11354577/files/output.tar.gz \
    && tar xzf output.tar.gz -C output \
    && rm output.tar.gz \
    && cd .. \
    && pip install --no-cache-dir ./abodybuilder3

# ---- Copy workflow node scripts ----------------------------------------
COPY . .

# Make shell scripts executable and normalize line endings
RUN find . -name "run.sh" -exec dos2unix {} \; \
    && find . -name "run.sh" -exec chmod +x {} \;
