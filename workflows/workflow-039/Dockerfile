# Antibody developability triage: ANARCI + ImmuneBuilder + AbLang2 + promb.
#
# This is a conda base image, not python:3.x-slim, and that is forced rather than preferred.
# ImmuneBuilder's refinement step needs OpenMM and pdbfixer, and its sequence checks
# `from anarci import validate_sequence, anarci, scheme_short_to_long` — classic ANARCI, which
# needs the HMMER binaries and the germline HMM library alongside the Python package. All three
# are conda-first (upstream's own install instructions say so), so this follows the
# workflow-016 pattern: conda for those, pip for everything else.
#
# CPU only, deliberately. Every node in this workflow runs on CPU: ABodyBuilder2 is
# ~a second per antibody and AbLang2 scores a panel in seconds, so a 500-clone run needs no
# GPU at all. That is the economic point of this workflow — see README "Requirements".
FROM continuumio/miniconda3

ENV PYTHONDONTWRITEBYTECODE=1 \
    TORCH_HOME=/models/torch \
    PROMB_CACHE=/models/promb \
    HF_HOME=/models/hf \
    OMP_NUM_THREADS=1

# Python 3.11: ANARCII (if ever swapped in, see the pipeline doc) requires >=3.11, promb >=3.8,
# and ImmuneBuilder is unpinned. 3.11 is the newest version all of conda-forge/bioconda carry
# these packages for.
RUN conda install -y -c conda-forge -c bioconda \
      python=3.11 \
      anarci \
      hmmer \
      openmm \
      pdbfixer \
    && conda clean -afy

# Torch CPU wheels only. The default index would pull the full CUDA runtime (~2.5 GiB) for a
# workflow that never touches a GPU.
RUN pip install --no-cache-dir \
      --index-url https://download.pytorch.org/whl/cpu \
      torch

# Pinned exactly, all BSD-3 or MIT:
#   ImmuneBuilder 1.2  — ABodyBuilder2 / NanoBodyBuilder2 / TCRBuilder2 (BSD-3)
#   ablang2      0.2.1 — paired antibody language model for the mutation scan (BSD-3)
#   promb        1.0.2 — human 9-mer peptide prevalence, the OASis-like score (MIT)
# biotite supplies SASA (Shrake-Rupley) for the hydrophobic/charged patch metrics, and
# biopython supplies the isoelectric point estimate. No LGPL/GPL components, and nothing
# academic-licence-only — see the pipeline doc on why NetMHCIIpan is absent.
RUN pip install --no-cache-dir \
      ImmuneBuilder==1.2 \
      ablang2==0.2.1 \
      promb==1.0.2 \
      "biotite>=1.5" \
      "biopython>=1.88" \
      numpy scipy pandas

# Bake every model weight and reference database into the image, so no node needs the network.
# ImmuneBuilder downloads its weights on first use (`are_weights_ready` / `download_file` in
# ImmuneBuilder/util.py), AbLang2 downloads on `pretrained()`, and promb fetches its predefined
# databases on `init_db`. Instantiating each one at build time is what populates the caches.
RUN python -c "\
from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2; \
ABodyBuilder2(); NanoBodyBuilder2(); \
print('ImmuneBuilder weights ready')"

RUN python -c "\
import ablang2; \
ablang2.pretrained(model_to_use='ablang2-paired', random_init=False, ncpu=1, device='cpu'); \
print('AbLang2 weights ready')"

RUN python -c "\
from promb import init_db; \
init_db('human-oas'); \
print('promb human-oas database ready')"

# Smoke-test ANARCI's HMM library at build time rather than discovering it is missing in node 01.
RUN python -c "\
from anarci import anarci; \
seq='EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS'; \
n,_,_=anarci([('h',seq)], scheme='imgt', output=False); \
assert n[0], 'ANARCI failed to number a known VH'; \
print('ANARCI ok')"

WORKDIR /workspace
