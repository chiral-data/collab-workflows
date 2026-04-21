FROM python:3.11

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    jq \
    gfortran \
    libopenblas-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir "numpy<2.0"
RUN CFLAGS="-Wno-incompatible-pointer-types" pip install --no-cache-dir lightdock

WORKDIR /workspace