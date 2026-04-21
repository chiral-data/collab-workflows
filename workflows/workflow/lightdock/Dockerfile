FROM python:3.11

RUN pip install lightdock

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace