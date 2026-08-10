FROM python:3.6.15-slim-bullseye

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Runtime and system dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Packaging tools compatible with Python 3.6.
RUN python -m pip install --no-cache-dir \
    pip==21.3.1 \
    setuptools==59.6.0 \
    wheel==0.37.1

# AMPlify runtime dependencies.
RUN python -m pip install --no-cache-dir \
    numpy==1.16.6 \
    scipy==1.5.4 \
    pandas==1.1.5 \
    scikit-learn==0.24.2 \
    biopython==1.79 \
    h5py==2.10.0 \
    absl-py==1.4.0 \
    gast==0.2.2 \
    grpcio==1.48.2 \
    protobuf==3.19.6 \
    tensorboard==1.12.2 \
    tensorflow==1.12.0 \
    keras==2.2.4

# Clone AMPlify into /opt and pin the upstream source version.
WORKDIR /opt

RUN git clone https://github.com/BirolLab/AMPlify.git /opt/AMPlify \
    && cd /opt/AMPlify \
    && git checkout 3a07713c25b8a21ef66d31d10e121989d26d9320 \
    && rm -rf /opt/AMPlify/.git

WORKDIR /workspace

# Smoke test used when the container is run without arguments.
CMD ["python", "/opt/AMPlify/src/AMPlify.py", "--help"]