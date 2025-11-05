# Base image with Miniconda
FROM continuumio/miniconda3:latest

# Create environment and install dependencies
RUN conda create -y -n prep_env python=3.10 \
    && conda install -y -n prep_env -c conda-forge pdbfixer biopython \
    && conda clean -afy

# Activate environment automatically
SHELL ["conda", "run", "-n", "prep_env", "/bin/bash", "-c"]

# Install MGLTools for prepare_receptor4.py
RUN apt-get update && apt-get install -y wget unzip \
    && wget https://ccsb.scripps.edu/download/262/ -O mgltools.zip || true \
    && unzip mgltools.zip -d /opt/mgltools || true \
    && chmod +x /opt/mgltools/bin/*

# Set working directory
WORKDIR /workspace

# Copy all job files into the container
COPY . /workspace

# Entry point script
ENTRYPOINT ["bash", "run.sh"]
