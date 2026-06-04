FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+=2.12.0+ds-3+b1 \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    "biopython>=1.81" \
    "primer3-py==2.3.0" \
    "requests>=2.31.0"

WORKDIR /workspace
