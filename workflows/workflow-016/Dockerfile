FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ gfortran && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    lightdock \
    numpy \
    scipy \
    prody \
    freesasa

WORKDIR /workspace
