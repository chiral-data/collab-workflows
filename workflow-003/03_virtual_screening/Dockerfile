# 03_virtual_screening/Dockerfile
FROM continuumio/miniconda3:latest

WORKDIR /work

# Create and prepare environment
RUN conda create -n vina_env python=3.10 -y && \
    conda install -n vina_env -y -c conda-forge -c bioconda \
      vina p2rank openbabel biopython pandas && \
    conda run -n vina_env pip install openpyxl && \
    conda clean -afy

SHELL ["conda", "run", "-n", "vina_env", "/bin/bash", "-lc"]

# Copy all scripts
COPY run.sh 03_virtual_screening.py rank_vina.py /work/
RUN chmod +x /work/run.sh

ENTRYPOINT ["bash", "run.sh"]