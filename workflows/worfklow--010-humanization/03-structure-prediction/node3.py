# ============================================================
# workflow-010-antibody-structure-prediction
# CPU-only image (Intel / AMD x86_64)
# Includes ProtT5 support for ABB3-LM (Node 02)
# ============================================================

FROM python:3.10-slim

# ---- System dependencies -----------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        curl \
        dos2unix \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workflow

# ---- Python dependencies -----------------------------------------------
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

# ProtT5 dependencies (transformers + sentencepiece for Node 02 / ABB3-LM)
RUN pip install --no-cache-dir \
    transformers \
    sentencepiece

# ---- ABodyBuilder3 -------------------------------------------------------
# Clone the upstream repo, download model checkpoints from Zenodo via
# the repo's own download.sh, then install as a Python package.
# This ensures checkpoints are always in sync with the cloned version.
RUN git clone https://github.com/oxpig/ABodyBuilder3.git abodybuilder3 \
    && cd abodybuilder3 \
    && bash download.sh \
    && cd .. \
    && pip install --no-cache-dir ./abodybuilder3

# ---- Copy workflow node scripts ----------------------------------------
COPY 01-fasta-validation/node1.py        ./01-fasta-validation/node1.py
COPY 01-fasta-validation/run.sh          ./01-fasta-validation/run.sh
COPY 02-plm-embedding/node2.py           ./02-plm-embedding/node2.py
COPY 02-plm-embedding/run.sh             ./02-plm-embedding/run.sh
COPY 03-structure-prediction/node3.py    ./03-structure-prediction/node3.py
COPY 03-structure-prediction/run.sh      ./03-structure-prediction/run.sh
COPY 04-visualization-report/node4.py   ./04-visualization-report/node4.py
COPY 04-visualization-report/run.sh     ./04-visualization-report/run.sh

# Make shell scripts executable and normalize line endings
RUN find . -name "run.sh" -exec dos2unix {} \; \
    && find . -name "run.sh" -exec chmod +x {} \;

# ---- Runtime defaults --------------------------------------------------
# All overridable at `docker run` time via -e flags.
ENV HEAVY_FASTA=""         \
    LIGHT_FASTA=""         \
    CHECKPOINT_PATH=""     \
    USE_PLM=0              \
    PRECOMPUTED_DIR=""     \
    DEVICE=cpu             \
    REPORT_TITLE="ABB3 Structure Predictions"

# Default: run Node 01 only (entrypoint for individual node execution)
CMD ["bash", "01-fasta-validation/run.sh"]