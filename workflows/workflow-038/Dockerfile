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
#     docker run --rm --gpus all bioemu2:2026_11_20 \
#       python -c "import torch; print(torch.cuda.get_arch_list())"
# If sm_70 is absent, the floor is 7.5 and the V100 claim has to come out of the README.

ENV HF_HOME=/models \
    BIOEMU2_MSA_CACHE=/models/colabfold_msa \
    PYTHONDONTWRITEBYTECODE=1 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

# gcc and libc6-dev are a **runtime** requirement, not a build one. Triton compiles its kernel
# launchers on first use and shells out to a C compiler to do it; without one, sampling fails
# inside EmuSession.sample() with "Failed to find C compiler. Please specify via CC environment
# variable". Measured on the image without them: 1UBQ (76 tokens), 5PTI (65 tokens, three
# non-polymer components) and 2ZTA (72 tokens) all sampled fine and never reached that path,
# while 1CLL (155 tokens) failed there after the model had loaded. So it is not about having
# ligands, and the first failure looked like a size limit or a bad input rather than a missing
# package. python:3.12-slim ships neither compiler nor libc headers.
RUN apt-get update && apt-get install -y --no-install-recommends \
      git curl ca-certificates gcc libc6-dev && \
    rm -rf /var/lib/apt/lists/*

# 0.1.0rc0 is a pre-release, so the version must be pinned exactly (or --pre passed) or pip
# will refuse to install it at all.
#
# openmm needs no conda base image, which an earlier TODO here assumed it might. bioemu2
# requires it unpinned, and PyPI publishes cp312 manylinux wheels for it — the image built on
# 2026-11-20 resolved openmm 8.6.1-cp312-cp312-manylinux_2_28_x86_64, and manylinux_2_28 needs
# glibc >= 2.28 while python:3.12-slim (Debian 12) has 2.36. `pip download openmm
# --python-version 312 --only-binary=:all: --platform manylinux_2_28_x86_64` reproduces it
# without building anything. Measured in the built image: openmm 8.6.1, mdtraj 1.11.2,
# torch 2.13.0+cu126, numpy 2.4.6, biotite 2.0.0, bioemu2 0.1.0rc0.
RUN pip install --no-cache-dir \
      --extra-index-url https://download.pytorch.org/whl/cu126 \
      "bioemu2[cuda12]==0.1.0rc0"

# Bake the 1.91 GiB bioemu2-preview checkpoint into the image at a pinned revision, so the
# sampling nodes run with HF_HUB_OFFLINE=1 and never reach for the network.
RUN python -c "from huggingface_hub import snapshot_download; \
    snapshot_download('microsoft/bioemu2', \
      revision='dd4b0c1b69984b68edcb6200a544213be4d75242', \
      allow_patterns=['checkpoints/bioemu2-preview/*'])"

# Pre-warm the chemical component definitions.
#
# Anything that is not a standard residue — an ion, a cofactor, a terminal cap, a modified
# residue — is built from its CCD entry, which bioemu2 downloads from files.rcsb.org the first
# time it meets the code (bioemu2/structure/templates.py: `_CCD_URL`) and caches, along with
# the template it derives from it. That download happens inside `system.topology`, which nodes
# 01, 03, 04 and 05 all call, so without this step the first run of any holo system reaches the
# network from four frames deep in the library and — with no egress — fails there with a DNS
# error that says nothing about chemistry. Measured on the fixture entries: 1CLL needs CA,
# 1A3N HEM, 1AKE AP5, 5PTI PO4, 2ZTA ACE; 1UBQ and 1BNA need none.
#
# The list is ions, caps, cofactors, crystallisation additives, sugars and the frequent
# modified residues: 126 codes, about 2 MB of cache, every one verified to build in this image.
# UNX, UNK and DUM are in the list because a depositor's placeholder for unidentified density
# is still a component BioEmu2 will sample; node 01 warns when one is present. UNL is the one
# code that cannot be modelled at all -- it raises inside bioemu2's own template builder -- so
# it is not here, and node 01 names the components and the exclude_ligands parameter instead.
# A code outside it still needs network on first use, and every node that builds a topology
# says exactly that when the fetch fails, naming the code -- which is how EOH, the ethanol in
# 1CLL, was found: the fixture suite runs with no network and the failure named the component.
ENV BIOEMU2_CCD_CACHE=/models/ccd \
    BIOEMU2_TEMPLATE_CACHE=/models/templates

RUN python -c "\
from bioemu2.structure.inference_inputs import inference_input_structure_from_ccd as build; \
codes = 'CA MG ZN NA K CL MN FE FE2 CU CO NI CD HG ACE NH2 FOR HOH HEM ATP ADP AMP GTP GDP \
AP5 NAD NAP NDP FAD FMN SAM SAH COA PLP SO4 PO4 GOL EDO PEG ACT CIT TRS MES EPE DMS IOD BR \
MSE SEP TPO PTR CSO CME KCX MLY HYP PCA ABA NLE ORN AIB \
EOH MPD PG4 P6G 1PE BME DTT FMT ACY NO3 CO3 NH4 AZI SCN CAC BCT TLA MLI OXL SIN \
NAG NDG MAN BMA GAL GLC FUC XYP SIA ANP ACP AGS GNP GSP SF4 FES HEC HEB CLA BCL F3S \
MOO WO4 VO4 ALF BEF PGE 2PE TBU IPA MOH ACN URE GSH BOG LDA HED IMD PYR \
XE KR UNX UNK DUM DOD'.split(); \
[(build(code), print('warmed', code, flush=True)) for code in codes]; \
print(f'pre-warmed {len(codes)} CCD codes')"

# Report rendering (node 06) and analysis (node 05) need nothing beyond BioEmu2's own deps:
# mdtraj, biotite, numpy and scipy all come in with bioemu2.

WORKDIR /workspace
