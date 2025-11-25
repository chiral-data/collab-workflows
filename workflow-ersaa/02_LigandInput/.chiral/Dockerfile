# Node 2: Download ligands from PubChem
FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

RUN conda create -y -n env python=3.10 \
    && conda install -y -n env -c conda-forge requests \
    && conda clean -afy

WORKDIR /workspace
COPY .. /workspace

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "env", "bash", "run.sh"]
