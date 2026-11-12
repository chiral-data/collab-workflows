FROM python:3.12-slim

# RFD4-Proteina is Linux x86_64 + Python 3.12 only. Upstream's own Dockerfile uses a slim
# Python base rather than a CUDA base image, selecting the accelerator at install time.
ENV HF_HOME=/models \
    DATA_PATH=/opt/rfd4/data_path \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

RUN apt-get update && apt-get install -y --no-install-recommends \
      git curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Upstream warns that a fresh venv's bundled pip trusts a compromised CA — pin pip first.
RUN pip install --no-cache-dir "pip==26.2.2" && pip install --no-cache-dir "uv>=0.13.4"

# CUDA 12 profile, matching this repo's GPU estate (no CUDA 13 host driver available).
# The BioIR backend is cu13-only and is deliberately forgone.
# TODO(build): pin the rfd4-proteina version once `uv pip index versions` is run at build time.
RUN uv pip install --system --no-cache --no-sources \
      --torch-backend cu126 "rfd4-proteina[all-cuda12]"

# Node dependencies: structure parsing (01, 04, 05), spec templating (02),
# .cif.zst inputs (01), RCSB fetch (01). Pinned so the CPU nodes are reproducible
# independently of whatever rfd4-proteina happens to pull in.
RUN uv pip install --system --no-cache \
      "biopython==1.85" \
      "jinja2==3.1.6" \
      "zstandard==0.23.0" \
      "requests==2.32.4"

# Bake the ~1.74 GB matched flow + autoencoder pair into the image at a pinned revision.
RUN python -c "from huggingface_hub import snapshot_download; \
    snapshot_download('nvidia/RFD4-Proteina', revision='a06e2509c7033b010fc43614e9fd840303588502')"

WORKDIR /workspace
