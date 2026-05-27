FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install Core & Scientific packages
RUN micromamba install -y -n base -c conda-forge \
    python=3.11 \
    numpy \
    pandas \
    openpyxl \
    matplotlib \
    scikit-learn \
    && micromamba clean -afy

# Install Cheminformatics tools (RDKit)
RUN micromamba install -y -n base -c conda-forge \
    rdkit \
    && micromamba clean -afy

# Install Deep Learning Framework (TensorFlow/Keras)
RUN micromamba install -y -n base -c conda-forge \
    tensorflow \
    keras \
    && micromamba clean -afy

# Fix TensorFlow execstack issue on kernels with noexecstack
USER root
RUN apt-get update -qq && apt-get install -y -q patchelf && \
    find /opt/conda -name "libtensorflow_cc.so*" -exec patchelf --clear-execstack {} \; && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

LABEL description="QSAR Workflow Environment"
LABEL tools="RDKit, TensorFlow, Scikit-learn"

CMD ["python", "-c", "import rdkit; import tensorflow; import sklearn; print('qsar ready')"]
