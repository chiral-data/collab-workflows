# Node1: download receptor
FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash","-c"]

# Install fast packages with mamba (as root)
RUN conda install -y -n base -c conda-forge mamba && \
    mamba install -y -n base -c conda-forge python=3.10 biopython requests && \
    conda clean -afy

# Switch to working directory BEFORE copy
WORKDIR /workspace

# Copy all files while still root
COPY . /workspace

# Give run.sh permission while still root
RUN chmod +x /workspace/run.sh

# Create non-root user and switch to it
RUN useradd -m -u 1000 worker
USER worker

# Final entrypoint
ENTRYPOINT ["bash", "/workspace/run.sh"]
