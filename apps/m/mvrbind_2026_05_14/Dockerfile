FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.10 \
    python3.10-dev \
    python3-pip \
    git \
    wget \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1

RUN pip install --upgrade pip

# PyTorch + CUDA 12.1
RUN pip install --no-cache-dir \
    torch==2.2.0 \
    torchvision==0.17.0 \
    torchaudio==2.2.0 \
    --index-url https://download.pytorch.org/whl/cu121

# torch-geometric sparse extensions
RUN pip install --no-cache-dir \
    torch-geometric==2.6.1 \
    torch-cluster==1.6.3 \
    torch-scatter==2.1.2 \
    torch-sparse==0.6.18 \
    torch-spline-conv==1.2.2 \
    -f https://data.pyg.org/whl/torch-2.2.0+cu121.html

# Python dependencies
RUN pip install --no-cache-dir \
    numpy==1.23.5 \
    biopython==1.85 \
    scikit-learn==1.7.0 \
    scipy==1.15.3 \
    matplotlib==3.10.3 \
    tqdm==4.67.1 \
    pyyaml==6.0.2

# Clone MVRBind (includes bundled ClustalW binary)
WORKDIR /opt
RUN git clone https://github.com/cschen-y/MVRBind.git

WORKDIR /opt/MVRBind
RUN chmod +x clustalw/clustalw2

WORKDIR /workspace

# Preprocessing pipeline (run in order before predict_binding_site.py):
#   Step 1: python /opt/MVRBind/data_process/get_pdb.py --pdb_id_chain <ID> --pdb_dir ./data/pdb
#   Step 2: python /opt/MVRBind/data_process/msa_score.py --clustalw /opt/MVRBind/clustalw/clustalw2 ...
#   Step 3: obtain secondary structure from http://rnapdbee.cs.put.poznan.pl/ (manual web step)
#   Step 4: generate RNABERT embeddings and RNAsnap2 SASA features (external tools)
CMD ["python", "/opt/MVRBind/predict_binding_site.py", "--help"]
