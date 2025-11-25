# Node 3: Prepare receptor using PDBFixer + MGLTools
FROM continuumio/miniconda3:24.1.2-0
SHELL ["/bin/bash", "-c"]

# Create conda environment and install dependencies
RUN conda create -y -n env python=3.10 \
    && conda install -y -n env -c conda-forge pdbfixer openmm biopython \
    && conda clean -afy

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends wget unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# MGLTools installation
RUN wget https://ccsb.scripps.edu/download/532/ -O /tmp/mgl.tgz \
    && mkdir -p /opt/mgltools \
    && tar -xzf /tmp/mgl.tgz -C /opt/mgltools --strip-components=1 \
    && rm /tmp/mgl.tgz

# Copy workspace and make run.sh executable
WORKDIR /workspace
COPY .. /workspace
RUN chmod +x run.sh

# Run the script inside the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "env", "bash", "run.sh"]
