FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget curl git build-essential \
    libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Copy environment and create conda env
COPY environment.yml /tmp/environment.yml
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda env create -f /tmp/environment.yml && conda clean -afy

# Make the workflow environment Python the default interpreter for non-interactive runners (e.g., Silva docker exec).
ENV CONDA_DEFAULT_ENV=diffdock_abag
ENV PATH=/opt/conda/envs/diffdock_abag/bin:/opt/conda/bin:$PATH

# Install DiffDock-PP and expose a stable path for workflow auto-detection.
RUN git clone --depth=1 https://github.com/ketatam/DiffDock-PP.git /opt/DiffDock-PP && \
    /opt/conda/envs/diffdock_abag/bin/pip install --no-cache-dir \
      dill \
      biopandas \
            "e3nn==0.5.1" \
      wandb \
      tensorboard \
            tensorboardX && \
        TORCH_MM=$(/opt/conda/envs/diffdock_abag/bin/python -c 'import torch; v=torch.__version__.split("+")[0].split("."); print(f"{v[0]}.{v[1]}.0")') && \
        CUDA_VER=$(/opt/conda/envs/diffdock_abag/bin/python -c 'import torch; print(torch.version.cuda or "cpu")') && \
        if [ "$CUDA_VER" = "cpu" ]; then CUDA_TAG=cpu; else CUDA_TAG=cu$(echo "$CUDA_VER" | tr -d '.'); fi && \
        PYG_URL="https://data.pyg.org/whl/torch-${TORCH_MM}+${CUDA_TAG}.html" && \
        /opt/conda/envs/diffdock_abag/bin/pip install --no-cache-dir \
            torch-scatter \
            torch-sparse \
            torch-cluster \
            torch-spline-conv \
            -f "$PYG_URL" && \
        /opt/conda/envs/diffdock_abag/bin/pip install --no-cache-dir torch-geometric
ENV DIFFDOCK_PP_PATH=/opt/DiffDock-PP

SHELL ["/bin/bash", "-c"]
RUN echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc && \
    echo 'conda activate diffdock_abag' >> ~/.bashrc

WORKDIR /workspace
CMD ["/bin/bash"]