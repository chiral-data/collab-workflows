# ============================================================
# ABB3 Antibody Structure Prediction Workflow
# CPU-only image (Intel / AMD)
# Includes ProtT5 support for ABB3-LM (Node 2)
# ============================================================

# Slim Python 3.10 base — no CUDA needed for CPU-only
FROM python:3.10-slim

# ---- System dependencies -----------------------------------------------
# build-essential: needed to compile some Python packages from source
# git: needed to install abodybuilder3 from GitHub
# wget: useful for downloading checkpoints/data
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        dos2unix \
    && rm -rf /var/lib/apt/lists/*

# ---- Working directory -------------------------------------------------
WORKDIR /workflow

# ---- Python dependencies -----------------------------------------------
# Install in a specific order to avoid conflicts.
# torch CPU-only wheel is much smaller than the default (no CUDA bundled).
RUN pip install --no-cache-dir --upgrade pip

# PyTorch CPU-only (avoids pulling in ~2 GB of CUDA libraries)
RUN pip install --no-cache-dir \
    torch==2.1.2 \
    --index-url https://download.pytorch.org/whl/cpu

# Core ML / bio dependencies
RUN pip install --no-cache-dir \
    numpy \
    scipy \
    biopython \
    ml-collections \
    einops

# ProtT5 dependencies (transformers + sentencepiece for Node 2 / ABB3-LM)
RUN pip install --no-cache-dir \
    transformers \
    sentencepiece

# abodybuilder3 — installed from a local clone of the repo
# (clone it first: git clone https://github.com/oxpig/ABodyBuilder3.git)
COPY abodybuilder3/ ./abodybuilder3/
RUN pip install --no-cache-dir ./abodybuilder3

# ---- Copy workflow scripts into the image ------------------------------
COPY node1/ ./node1/
COPY node2/ ./node2/
COPY node3/ ./node3/
COPY node4/ ./node4/
COPY run1.sh run2.sh run3.sh run4.sh run_all.sh ./

# Make shell scripts executable and fix Windows line endings
RUN dos2unix run1.sh run2.sh run3.sh run4.sh run_all.sh && \
    chmod +x run1.sh run2.sh run3.sh run4.sh run_all.sh

# ---- Runtime defaults --------------------------------------------------
# These can all be overridden at `docker run` time via -e flags.
ENV HEAVY_FASTA=data/heavy.fasta \
    LIGHT_FASTA=data/light.fasta \
    CHECKPOINT_PATH="" \
    USE_PLM=0 \
    PRECOMPUTED_DIR="" \
    DEVICE=cpu \
    REPORT_TITLE="ABB3 Structure Predictions" \
    RESULTS_ROOT=results

# Default command — runs the full pipeline.
# Override with e.g. `docker run ... bash run3.sh` to run a single node.
CMD ["bash", "run_all.sh"]