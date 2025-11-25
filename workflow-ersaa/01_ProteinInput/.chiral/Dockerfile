# Node1: download receptor
FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash","-c"]

# install fast packages with mamba
RUN conda install -y -n base -c conda-forge mamba && \
    mamba install -y -n base -c conda-forge python=3.10 biopython requests && \
    conda clean -afy

RUN useradd -m -u 1000 worker
USER worker
WORKDIR /workspace

COPY . /workspace
RUN chmod +x /workspace/run.sh

ENTRYPOINT ["bash","/workspace/run.sh"]
