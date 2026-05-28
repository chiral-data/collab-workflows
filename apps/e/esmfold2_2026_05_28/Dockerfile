FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Python 3.12 via deadsnakes PPA (esm requires >=3.12,<3.13)
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common ca-certificates gnupg \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-dev python3.12-venv \
    git build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN python3.12 -m ensurepip --upgrade \
    && python3.12 -m pip install --upgrade pip setuptools wheel \
    && update-alternatives --install /usr/bin/python  python  /usr/bin/python3.12 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1 \
    && update-alternatives --install /usr/bin/pip     pip     /usr/local/bin/pip3.12 1

# PyTorch 2.3 + CUDA 12.1 (install before esm to ensure CUDA build is used)
RUN pip install --no-cache-dir \
    torch==2.3.0 \
    --index-url https://download.pytorch.org/whl/cu121

# ESMFold2 package — pulls in Biohub transformers fork + all declared deps
RUN pip install --no-cache-dir \
    "esm @ git+https://github.com/Biohub/esm.git@c94ed8d"

# Pre-download model weights at build time (token passed via --secret, never baked in)
ENV HF_HOME=/models
RUN --mount=type=secret,id=hf_token \
    huggingface-cli download biohub/ESMFold2 \
      --token "$(cat /run/secrets/hf_token)" --local-dir-use-symlinks False && \
    huggingface-cli download Biohub/ESMC-6B \
      --token "$(cat /run/secrets/hf_token)" --local-dir-use-symlinks False

# Smoke test (requires --gpus all):
#   docker run --gpus all <image> python -c "
#     from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
#     model = ESMFold2Model.from_pretrained('biohub/ESMFold2').cuda().eval()
#     print('ESMFold2 ready')
#   "

WORKDIR /workspace

CMD ["python", "-c", "from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model; print('ESMFold2 ready')"]
