# Node 3: Prepare receptor using PDBFixer + MGLTools
FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash", "-c"]

# Conda env for preparation
RUN conda create -y -n prep_env python=3.10 \
    && conda install -y -n prep_env -c conda-forge pdbfixer openmm biopython \
    && conda clean -afy

# System deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget unzip && apt-get clean && rm -rf /var/lib/apt/lists/*

# MGLTools
RUN wget https://ccsb.scripps.edu/download/532/ -O /tmp/mgltools.tar.gz \
    && mkdir -p /opt/mgltools \
    && tar -xzf /tmp/mgltools.tar.gz -C /opt/mgltools --strip-components=1 \
    && rm /tmp/mgltools.tar.gz

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "prep_env", "bash", "run.sh"]
