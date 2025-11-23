# Node 1: Download PDB receptor
FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

# Create environment
RUN conda create -y -n dl_env python=3.10 \
    && conda install -y -n dl_env -c conda-forge biopython requests \
    && conda clean -afy

# Create non-root user
RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x /workspace/run.sh

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "dl_env", "bash", "run.sh"]
