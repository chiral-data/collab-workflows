FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash","-c"]

# Use system Java and python libs
RUN apt-get update && apt-get install -y --no-install-recommends openjdk-11-jre-headless wget && \
    conda install -y -n base -c conda-forge mamba && \
    mamba install -y -n base -c conda-forge python=3.10 pandas && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && conda clean -afy

RUN useradd -m -u 1000 worker
USER worker
WORKDIR /workspace

COPY . /workspace
RUN chmod +x /workspace/run.sh

ENTRYPOINT ["bash","/workspace/run.sh"]
