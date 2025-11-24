# Node 4: Prepare ligands using MGLTools + OpenBabel
FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash", "-c"]

RUN conda create -y -n ligprep_env python=3.10 \
    && conda install -y -n ligprep_env -c conda-forge openbabel \
    && conda clean -afy

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

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ligprep_env", "bash", "run.sh"
