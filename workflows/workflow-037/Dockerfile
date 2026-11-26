FROM python:3.12-slim

# RFD4-Proteina is Linux x86_64 + Python 3.12 only. Upstream's own Dockerfile uses a slim
# Python base rather than a CUDA base image, selecting the accelerator at install time.
ENV HF_HOME=/models \
    DATA_PATH=/opt/rfd4/data_path \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

# libxrender1/libxext6/libsm6 are not optional. Without them Open Babel's plugin load
# fails, GetSupportedInputFormat() returns malformed data, and atomworks.ml's import of it
# dies with "not enough values to unpack" — which takes rfproteina's condition_spec, and
# therefore node 02's validation and generation itself, down with it.
# gcc/libc6-dev are needed at *run* time, not build time: Triton compiles its kernels on
# the fly and aborts with "Failed to find C compiler" without one. This happens even with
# compile_lengths=null, which only skips the torch.compile warm-up, so a compiler is not
# optional for generation. Upstream's own Dockerfile does not mention it because its base
# image already carries one.
RUN apt-get update && apt-get install -y --no-install-recommends \
      git curl ca-certificates gcc libc6-dev libxrender1 libxext6 libsm6 && \
    rm -rf /var/lib/apt/lists/*

# Upstream warns that a fresh venv's bundled pip trusts a compromised CA — pin pip first.
RUN pip install --no-cache-dir "pip==26.2.2" && pip install --no-cache-dir "uv>=0.13.4"

# CUDA 12 profile, matching this repo's GPU estate (no CUDA 13 host driver available).
# The BioIR backend is cu13-only and is deliberately forgone.
# Pinned: 0.0.3 is the newest release (checked 2026-11-12); requires-python is >=3.12,<3.13.
RUN uv pip install --system --no-cache --no-sources \
      --torch-backend cu126 "rfd4-proteina[all-cuda12]==0.0.3"

# Node dependencies: structure parsing (01, 04, 05), spec emission (02),
# .cif.zst inputs (01), RCSB fetch (01). Pinned so these do not drift with whatever
# rfd4-proteina happens to pull in. Not a complete pin of the CPU nodes: node 04
# imports numpy, which is not pinned here and comes from rfd4-proteina's closure.
RUN uv pip install --system --no-cache \
      "biopython==1.85" \
      "pyyaml==6.0.2" \
      "zstandard==0.23.0" \
      "requests==2.32.4"

# Bake the 1.74 GB model bundle in at a pinned revision.
#
# The patterns are not an optimisation, they are what the loader asks for:
# rfproteina/utils/hub.py:load_pretrained_bundle() calls snapshot_download with
# allow_patterns=["config.json", "model.safetensors"], and the autoencoder config travels
# inside that bundle. The repository as a whole is 82.57 GB, because it also ships the
# AlphaFold2 (51.6 GB) and OpenFold3 (21.2 GB) refolding runtime images used for scoring,
# plus native .ckpt copies of the weights. This workflow runs ~metrics and never scores,
# so none of that is fetched. An unfiltered snapshot_download fills the disk and fails.
#
# snapshot_download writes its tree cache 0600, so a non-root reader cannot open it,
# huggingface_hub reports it corrupted and falls back to the network -- which
# HF_HUB_OFFLINE=1 then refuses, failing generation with a misleading offline-mode error.
# Read and traverse only: nothing here is written at run time. The chmod shares this RUN
# because in a later one it would copy all 1.74 GB of weights into a new layer.
RUN python -c "from huggingface_hub import snapshot_download; \
    snapshot_download('nvidia/RFD4-Proteina', \
    revision='a06e2509c7033b010fc43614e9fd840303588502', \
    allow_patterns=['config.json', 'model.safetensors'])" && \
    chmod -R a+rX /models

# Everything above is baked as root, but silva runs containers as the invoking host's
# UID to keep bind-mounted outputs owned by the user. That UID has no /etc/passwd entry,
# so HOME resolves to "/" and every library that caches under $HOME dies on a read-only
# root: Triton's kernel cache takes generation down with it, matplotlib's config dir only
# warns. Point both at a writable path.
ENV HOME=/tmp \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/.config/matplotlib

WORKDIR /workspace
