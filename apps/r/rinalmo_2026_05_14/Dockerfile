FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3-pip \
    git \
    wget \
    ninja-build \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

RUN pip install --upgrade pip setuptools wheel packaging

# PyTorch 2.1.0 + CUDA 11.8 (matches environment.yml); numpy pinned to 1.x
# because PyTorch 2.1.0 was compiled against NumPy 1.x ABI
RUN pip install --no-cache-dir \
    torch==2.1.0 \
    torchvision==0.16.0 \
    torchaudio==2.1.0 \
    numpy==1.26.0 \
    --index-url https://download.pytorch.org/whl/cu118

# flash-attn 2.3.2 — requires CUDA devel headers + ninja; packaging must be
# pre-installed because --no-build-isolation skips the isolated build env
RUN pip install --no-cache-dir flash-attn==2.3.2 --no-build-isolation

# Clone and install RiNALMo (installs: torch, numpy, ml_collections, packaging, gdown)
WORKDIR /opt
RUN git clone https://github.com/lbcb-sci/RiNALMo.git

WORKDIR /opt/RiNALMo
RUN pip install --no-cache-dir .

# Additional runtime dependencies from environment.yml
RUN pip install --no-cache-dir \
    einops==0.7.0 \
    biopython==1.81 \
    pytorch-lightning==2.1.3 \
    torchmetrics==1.2.1

# Model weights (~2.5 GB for Giga) are NOT baked in.
# Download at runtime via gdown, or mount a weights volume:
#   docker run --gpus all -v /path/to/weights:/weights ...
# Weights are auto-downloaded on first get_pretrained_model() call to ~/.cache/rinalmo/

WORKDIR /workspace

# Usage:
#   from rinalmo.pretrained import get_pretrained_model
#   model, alphabet = get_pretrained_model(model_name="giga-v1")   # 650M
#   model, alphabet = get_pretrained_model(model_name="mega-v1")   # 150M
#   model, alphabet = get_pretrained_model(model_name="micro-v1")  # 35M
CMD ["python", "-c", "from rinalmo.pretrained import get_pretrained_model; print('RiNALMo ready')"]
