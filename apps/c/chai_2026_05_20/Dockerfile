# =============================================================================
# Chai-Lab Docker Image
# Multi-stage build with pre-baked model weights
# Based on upstream Dockerfile.chailab
# =============================================================================

ARG CHAI_VERSION=0.6.1

# ---------------------------------------------------------------------------
# Stage 1: Builder — install chai_lab and download model weights
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS builder

ENV LD_LIBRARY_PATH=""
ENV \
  LANG=C.UTF-8 \
  LC_ALL=C.UTF-8 \
  DEBIAN_FRONTEND=noninteractive \
  PYTHONUNBUFFERED=TRUE \
  PYTHONFAULTHANDLER=1 \
  PYTHONPYCACHEPREFIX='/tmp/.chai_lab_pycache' \
  MYPY_CACHE_DIR='/tmp/.chai_lab_mypy_cache' \
  LD_LIBRARY_PATH="/usr/local/cuda/lib64/stubs/:${LD_LIBRARY_PATH}" \
  UV_HTTP_TIMEOUT=1000 \
  VIRTUAL_ENV=/opt/venv \
  CHAI_DOWNLOADS_DIR=/opt/chai-weights

ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN --mount=type=cache,target=/var/cache/apt \
  apt-get -qq update \
  && apt-get -qq install -y \
  gnupg ca-certificates wget git curl \
  clang build-essential libstdc++6 \
  python3.10 python3.10-dev \
  kalign \
  && ln -s /lib/x86_64-linux-gnu/libcuda.so.1 /lib/x86_64-linux-gnu/libcuda.so || true \
  && ldconfig /lib/x86_64-linux-gnu/

ARG CHAI_VERSION
RUN --mount=type=cache,target=/root/.cache/uv \
  curl -LsSf https://astral.sh/uv/0.5.4/install.sh | sh \
  && . $HOME/.local/bin/env \
  && uv venv --no-python-downloads $VIRTUAL_ENV \
  && . $VIRTUAL_ENV/bin/activate \
  && uv pip install uv pip chai_lab==${CHAI_VERSION}

# Pre-download all model weights and data via curl (no GPU needed)
ARG CHAI_ASSETS_URL=https://chaiassets.com/chai1-inference-depencencies
RUN mkdir -p $CHAI_DOWNLOADS_DIR/models_v2 $CHAI_DOWNLOADS_DIR/esm \
  && echo "Downloading conformers..." \
  && curl -fSL -o $CHAI_DOWNLOADS_DIR/conformers_v1.apkl \
     ${CHAI_ASSETS_URL}/conformers_v1.apkl \
  && echo "Downloading ESM2 model..." \
  && curl -fSL -o $CHAI_DOWNLOADS_DIR/esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt \
     ${CHAI_ASSETS_URL}/esm2/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt \
  && echo "Downloading model components..." \
  && for comp in feature_embedding.pt bond_loss_input_proj.pt token_embedder.pt trunk.pt diffusion_module.pt confidence_head.pt; do \
       echo "  -> $comp" \
       && curl -fSL -o $CHAI_DOWNLOADS_DIR/models_v2/$comp \
          ${CHAI_ASSETS_URL}/models_v2/$comp; \
     done \
  && echo "All weights downloaded successfully"

# ---------------------------------------------------------------------------
# Stage 2: Runtime — slim image with only what's needed
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS runtime

ARG CHAI_VERSION
LABEL org.opencontainers.image.title="chai-lab" \
      org.opencontainers.image.version="${CHAI_VERSION}" \
      org.opencontainers.image.description="Chai-1 molecular structure prediction"

ENV LD_LIBRARY_PATH=""
ENV \
  LANG=C.UTF-8 \
  LC_ALL=C.UTF-8 \
  DEBIAN_FRONTEND=noninteractive \
  PYTHONUNBUFFERED=TRUE \
  PYTHONFAULTHANDLER=1 \
  PYTHONPYCACHEPREFIX='/tmp/.chai_lab_pycache' \
  LD_LIBRARY_PATH="/usr/local/cuda/lib64/stubs/:${LD_LIBRARY_PATH}" \
  VIRTUAL_ENV=/opt/venv \
  CHAI_DOWNLOADS_DIR=/opt/chai-weights

ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Runtime-only system deps (no build-essential, clang, dev headers)
RUN apt-get -qq update \
  && apt-get -qq install -y --no-install-recommends \
  ca-certificates \
  libstdc++6 \
  python3.10 \
  kalign \
  libibverbs1 librdmacm1 \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && ln -s /usr/bin/python3.10 /usr/bin/python3 \
  && ln -s /lib/x86_64-linux-gnu/libcuda.so.1 /lib/x86_64-linux-gnu/libcuda.so || true \
  && ldconfig /lib/x86_64-linux-gnu/

# Copy virtualenv with all installed packages
COPY --from=builder /opt/venv /opt/venv

# Copy pre-downloaded model weights
COPY --from=builder /opt/chai-weights /opt/chai-weights

# Persist env vars for all shells
RUN echo "PATH=\"$PATH\"" >> /etc/environment \
  && echo "LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH\"" >> /etc/environment

WORKDIR /workspace

ENTRYPOINT []
CMD ["tail", "-f", "/dev/null"]
