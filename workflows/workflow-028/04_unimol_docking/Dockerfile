FROM dptechnology/unicore:latest-pytorch1.12.1-cuda11.6-rdma

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir rdkit biopandas

# Clone Uni-Mol repository for the docking V2 predict script and data pipeline
RUN git clone --depth 1 https://github.com/deepmodeling/Uni-Mol /opt/unimol \
    && pip install --no-cache-dir -e /opt/unimol/unimol_tools

RUN mkdir -p /opt/unimol_weights \
    && wget -q -O /opt/unimol_weights/unimol_docking_v2_240517.pt \
       "https://www.dropbox.com/scl/fi/sfhrtx1tjprce18wbvmdr/unimol_docking_v2_240517.pt?rlkey=5zg7bh150kcinalrqdhzmyyoo&st=n6j0nt6c&dl=1"

WORKDIR /workspace
