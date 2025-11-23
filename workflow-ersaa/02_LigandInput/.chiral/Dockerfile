# Node 2: Download ligands from PubChem
FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

RUN conda create -y -n dl_env python=3.10 \
    && conda install -y -n dl_env -c conda-forge requests \
    && conda clean -afy

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "dl_env", "bash", "run.sh"]
