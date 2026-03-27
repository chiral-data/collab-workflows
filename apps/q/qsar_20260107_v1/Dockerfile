# QSAR Workflow Environment
# Base image for quantitative structure-activity relationship modeling

FROM continuumio/miniconda3:latest

# Update conda and install mamba
RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

# Install Core & Scientific packages
RUN mamba install -y -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    openpyxl \
    matplotlib \
    scikit-learn \
    && mamba clean -afy

# Install Cheminformatics tools (RDKit)
RUN mamba install -y -c conda-forge \
    rdkit \
    && mamba clean -afy

# Install Deep Learning Framework (TensorFlow/Keras)
RUN mamba install -y -c conda-forge \
    tensorflow \
    keras \
    && mamba clean -afy

WORKDIR /workspace

# Metadata
LABEL description="QSAR Workflow Environment"
LABEL tools="RDKit, TensorFlow, Scikit-learn"

CMD ["python", "-c", "import rdkit; import tensorflow; import sklearn; print('qsar ready')"]
