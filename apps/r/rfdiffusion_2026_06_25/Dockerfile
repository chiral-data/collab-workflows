# RFdiffusion (original) — de novo protein binder design
# Based on official Dockerfile: github.com/RosettaCommons/RFdiffusion/blob/main/docker/Dockerfile
FROM nvcr.io/nvidia/cuda:11.6.2-cudnn8-runtime-ubuntu20.04

RUN apt-get -q update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
       git python3.9 python3-pip wget \
    && python3.9 -m pip install -q -U --no-cache-dir pip \
    && rm -rf /var/lib/apt/lists/*

# Clone RFdiffusion repository
RUN git clone https://github.com/RosettaCommons/RFdiffusion.git /app/RFdiffusion

WORKDIR /app/RFdiffusion

# Install Python dependencies (pinned versions from official Dockerfile)
RUN pip install -q --no-cache-dir \
    dgl==1.0.2+cu116 -f https://data.dgl.ai/wheels/cu116/repo.html \
    torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116 \
    "numpy<2" \
    e3nn==0.3.3 wandb==0.12.0 pynvml==11.0.0 \
    git+https://github.com/NVIDIA/dllogger#egg=dllogger \
    decorator==5.1.0 hydra-core==1.3.2 pyrsistent==0.19.3 \
    /app/RFdiffusion/env/SE3Transformer \
    && pip install --no-cache-dir /app/RFdiffusion --no-deps

# Download model weights for binder design
RUN mkdir -p models && \
    wget -q http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt -O models/Base_ckpt.pt && \
    wget -q http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt -O models/Complex_base_ckpt.pt && \
    wget -q http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt -O models/Complex_Fold_base_ckpt.pt

ENV DGLBACKEND="pytorch"

CMD ["python3.9", "scripts/run_inference.py", "--help"]
