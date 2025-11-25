# Node 1: Download PDB receptor
FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

RUN conda create -y -n env python=3.10 \
    && conda install -y -n env -c conda-forge biopython requests \
    && conda clean -afy

WORKDIR /workspace
COPY .. /workspace

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "env", "bash", "run.sh"]
