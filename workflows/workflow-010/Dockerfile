# ============================================================
# workflow-010-antibody-structure-prediction
# CUDA 12.1 image with ProtT5 support for ABB3-LM
# ============================================================

FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# ---- Environment -----------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV HF_HOME=/root/.cache/huggingface
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
    "transformers<=4.51.3" \
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

# ---- Pre-download ProtT5 weights (used by Node 02 when USE_PLM=1) ------
# Half-precision encoder-only variant: 2.4 GB instead of 11 GB
# Pass --build-arg HF_TOKEN=<token> to avoid HuggingFace anonymous rate-limiting
ARG HF_TOKEN=""
RUN HF_TOKEN=${HF_TOKEN} python3 -c "\
from transformers import T5EncoderModel, T5Tokenizer; \
T5EncoderModel.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc'); \
T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc')" \
    && chmod 755 /root \
    && chmod -R a+rX /root/.cache/huggingface \
    && rm -f /root/.cache/huggingface/token

# Use cached weights only — no network calls at runtime
ENV HF_HUB_OFFLINE=1

# ---- Copy workflow node scripts ----------------------------------------
COPY . .

# Make shell scripts executable and normalize line endings
RUN find . -name "run.sh" -exec dos2unix {} \; \
    && find . -name "run.sh" -exec chmod +x {} \;
