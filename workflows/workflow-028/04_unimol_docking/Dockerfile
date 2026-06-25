FROM dptechnology/unicore:latest-pytorch1.12.1-cuda11.6-rdma

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir rdkit biopandas

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
