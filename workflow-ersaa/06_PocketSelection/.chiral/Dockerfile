# Node 6: Select pocket and generate grids
FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

RUN conda create -y -n grid_env python=3.10 \
    && conda install -y -n grid_env -c conda-forge pandas \
    && conda clean -afy

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "grid_env", "bash", "run.sh"]

ENTRYPOINT ["bash", "run.sh"]
