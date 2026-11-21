FROM python:3.12-slim

# BioEmu2 is Linux x86_64 + NVIDIA GPU only, and requires Python >= 3.12. Upstream evaluated
# its benchmark metrics on 3.12 specifically, so this image stays on 3.12 rather than 3.13/3.14.
#
# The CUDA 12.6 install path is deliberate, for two reasons taken from bioemu2's own metadata.
#
# `[cuda12]` requires torch==2.13.0+cu126 exactly, plus cuequivariance{,-torch}==0.11.1 and the
# cu12 ops. `[cuda]` requires the same cuEquivariance version with the **cu13** ops and pins no
# torch at all, so torch resolves anywhere in bioemu2's own `>=2.11,<2.14` range and the image
# is a different image every time it is built. A pinned torch is the difference between a
# reproducible sampling run and one that cannot be reproduced after a point release.
#
# Second, cu126 runs on older drivers than the cu13 ops require, so this image is portable
# across a mixed estate rather than only to the newest host.
#
# An earlier version of this comment said the estate "tops out at 12.6 (see
# apps/applications.json)". That was wrong twice over: `cuda_version` in applications.json
# describes *images*, not host drivers, and the GPU test host's driver is 580.178.04, which
# supports CUDA 13. The reason to prefer cu126 is reproducibility and portability, not a
# missing driver.
#
# The cost is the CUDA 13 cuEquivariance kernels, which are a speed optimisation, not a
# capability. VERIFY ON THE GPU HOST, because the README's compute-capability floor depends on
# it and it has not been measured:
#     docker run --rm --gpus all bioemu2:2026_11_11 \
#       python -c "import torch; print(torch.cuda.get_arch_list())"
# If sm_70 is absent, the floor is 7.5 and the V100 claim has to come out of the README.

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
#
# openmm needs no conda base image, which an earlier TODO here assumed it might. bioemu2
# requires it unpinned, and PyPI publishes cp312 manylinux wheels for it — resolved 2026-11-20
# as openmm 8.3.1-cp312-cp312-manylinux_2_28_x86_64, and manylinux_2_28 needs glibc >= 2.28
# while python:3.12-slim (Debian 12) has 2.36. `pip download openmm --python-version 312
# --only-binary=:all: --platform manylinux_2_28_x86_64` reproduces it without building anything.
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
