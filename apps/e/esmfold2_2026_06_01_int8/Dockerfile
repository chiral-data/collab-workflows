FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3-pip \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

RUN pip install --upgrade pip setuptools wheel

# PyTorch 2.3 + CUDA 12.1 (install before esm to ensure CUDA build is used)
RUN pip install --no-cache-dir \
    torch==2.3.0 \
    --index-url https://download.pytorch.org/whl/cu121

# ESMFold2 package — pulls in Biohub transformers fork + all declared deps
RUN pip install --no-cache-dir \
    "esm@git+https://github.com/Biohub/esm.git@c94ed8d"

# TorchAO for int8 weight-only quantization of ESMC-6B backbone (~12 GB → ~8.5 GB)
RUN pip install --no-cache-dir torchao

# TODO: add weight pre-download once local inference is verified on GPU server
# Smoke test (requires ≥16 GB VRAM):
#   docker run --gpus all <image> python -c "
#     import torch
#     from torchao.quantization import Int8WeightOnlyConfig, quantize_
#     from transformers import ESMFold2Model
#     model = ESMFold2Model.from_pretrained('biohub/ESMFold2-Fast').eval()
#     quantize_(model, Int8WeightOnlyConfig())
#     model = model.cuda()
#     print('Peak VRAM:', torch.cuda.max_memory_allocated() / 1e9, 'GB')
#   "

WORKDIR /workspace

CMD ["python", "-c", "import esm; import torchao; print('esmfold2-int8 ready')"]
