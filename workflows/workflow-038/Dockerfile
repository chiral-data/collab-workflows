FROM python:3.12-slim

# BioEmu2 is Linux x86_64 + NVIDIA GPU only, and requires Python >= 3.12. Upstream evaluated
# its benchmark metrics on 3.12 specifically, so this image stays on 3.12 rather than 3.13/3.14.
#
# The CUDA 12.6 install path is deliberate. The `[cuda]` extra pulls cuequivariance-ops-cu13,
# i.e. a CUDA 13 driver, and this repo's GPU estate tops out at 12.6 (see apps/applications.json).
# `[cuda12]` pins torch==2.13.0+cu126 and the cu12 cuEquivariance ops, and is also the only path
# that supports compute capability 7.0 — so this image runs on V100 as well as T4/A100/L40S/H100.
# The cost is the CUDA 13 cuEquivariance kernels, which are a speed optimisation, not a capability.

ENV HF_HOME=/models \
    BIOEMU2_MSA_CACHE=/models/colabfold_msa \
    PYTHONDONTWRITEBYTECODE=1 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

RUN apt-get update && apt-get install -y --no-install-recommends \
      git curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# 0.1.0rc0 is a pre-release, so the version must be pinned exactly (or --pre passed) or pip
# will refuse to install it at all.
# TODO(build): openmm is a hard BioEmu2 dependency and is historically conda-first. If pip
# cannot resolve a cp312 wheel for it here, switch this base image to continuumio/miniconda3,
# `conda install -c conda-forge openmm` first, then pip-install bioemu2 (the workflow-016 pattern).
RUN pip install --no-cache-dir \
      --extra-index-url https://download.pytorch.org/whl/cu126 \
      "bioemu2[cuda12]==0.1.0rc0"

# Bake the 1.91 GiB bioemu2-preview checkpoint into the image at a pinned revision, so the
# sampling nodes run with HF_HUB_OFFLINE=1 and never reach for the network.
RUN python -c "from huggingface_hub import snapshot_download; \
    snapshot_download('microsoft/bioemu2', \
      revision='dd4b0c1b69984b68edcb6200a544213be4d75242', \
      allow_patterns=['checkpoints/bioemu2-preview/*'])"

# Report rendering (node 06) and analysis (node 05) need nothing beyond BioEmu2's own deps:
# mdtraj, biotite, numpy and scipy all come in with bioemu2.

WORKDIR /workspace
