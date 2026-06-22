FROM dptechnology/unicore:0.0.1-pytorch1.11.0-cuda11.3

# Uni-Core custom CUDA C++ extensions are pre-compiled in this base image
# for CUDA 11.3 / PyTorch 1.11.0. Rebuilding on a different CUDA version
# requires recompiling Uni-Core from source.

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir unimol_tools rdkit

# Clone Uni-Mol repository for the docking V2 predict script and data pipeline
RUN git clone --depth 1 https://github.com/deepmodeling/Uni-Mol /opt/unimol \
    && pip install --no-cache-dir -e /opt/unimol/unimol_tools

# ── Model weights ──────────────────────────────────────────────────────────
# Uni-Mol Docking V2 weights (464 MB) are NOT distributable via Docker.
# No pip or HuggingFace distribution exists.
#
# Before running this container:
#   1. Download the weights from the link in the Uni-Mol Docking V2 README:
#      https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_docking_v2
#   2. Place the downloaded files at /opt/unimol_weights/ on the host.
#   3. Set PARAM_WEIGHTS_PATH=/opt/unimol_weights in the job params, or
#      mount the directory into the container at that path.
# ──────────────────────────────────────────────────────────────────────────

WORKDIR /workspace
