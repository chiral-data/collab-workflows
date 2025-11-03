# 01_download/Dockerfile
FROM python:3.10-slim

LABEL maintainer="you@example.com"
ENV DEBIAN_FRONTEND=noninteractive

# Install system deps we'll need (wget, git, ca-certificates)
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
       wget \
       ca-certificates \
       git \
       bash \
 && rm -rf /var/lib/apt/lists/*

# Python deps
RUN pip install --no-cache-dir biopython requests tqdm

# Create work dir and set it as default
WORKDIR /workspace

# Copy all workflow files into container
COPY . /workspace

# Make shell scripts and python entry executable (defensive)
# RUN chmod +x /workspace/run.sh || true \
#  && chmod +x /workspace/download_job.py || true

# Set a default command (exec form) using python - more robust than relying on shebang
# ENTRYPOINT ["python", "/workspace/download_job.py"]
CMD ["/bin/bash"]
