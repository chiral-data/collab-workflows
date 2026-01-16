# Dockerfile for QSAR Workflow
# Inspired by Autodock Vina workflow architecture

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

# Copy workflow files
COPY . /workspace

# Ensure scripts are executable
RUN chmod +x /workspace/run_workflow.sh \
    /workspace/01_Data_Preparation/run_data_prep.sh \
    /workspace/02_Feature_Engineering/run_feature_eng.sh \
    /workspace/03_Model_Training/run_training.sh \
    /workspace/04_Prediction/run_prediction.sh

# Default command
CMD ["bash", "run_workflow.sh"]
