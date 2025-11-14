# Dockerfile for smina (Scoring and Minimization with AutoDock Vina)
# https://sourceforge.net/projects/smina/
#
# 2025-11-05: add mgltools (by Qin Wan)
# 2025-10-17: Modified by Qin Wan (rogerwq@gmail.com) 
# Created by Koki Shinbara 
#
# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    TZ=Etc/UTC

# システムパッケージのインストール
RUN apt-get update && apt-get install -y --no-install-recommends \
    # ビルド基盤
    build-essential wget curl git cmake g++ gcc pkg-config autoconf automake libtool \
    # 数値/科学計算
    liblapack-dev libblas-dev libgfortran5 libfftw3-dev libgsl-dev \
    libhdf5-dev libnetcdf-dev \
    # 圧縮/システム
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses-dev libreadline-dev libsqlite3-dev libdb-dev libgdbm-dev \
    libc6-dev libffi-dev libssl-dev libnss3-dev \
    # XML/画像まわり（描画はしないが依存回避用に最低限）
    libxml2-dev libxslt1-dev \
    libfreetype6-dev libpng-dev libjpeg-dev libxft-dev libxrender-dev libfontconfig1-dev \
    # OpenGL ランタイム（ヘッドレスでも一部ツールが要求）
    libgl1 libglu1-mesa \
    # 並列計算（MPI）
    libopenmpi-dev openmpi-bin \
    curl \ 
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# conda環境の設定
ENV CONDA_DEFAULT_ENV=base
ENV PATH="/opt/conda/bin:${PATH}"

# conda環境の作成とパッケージのインストール
RUN conda create -n in_silico python=3.13 -y

# 環境をアクティベート
ENV CONDA_DEFAULT_ENV=in_silico
ENV PATH="/opt/conda/envs/in_silico/bin:${PATH}"

# conda-forgeチャンネルを追加
RUN conda config --add channels conda-forge

# 主要なパッケージをcondaでインストール（段階的に）
# 基本的な科学計算パッケージ
RUN conda install -n in_silico -y \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    scikit-image

# 生物学・化学パッケージ
RUN conda install -n in_silico -y \
    biopython \
    openmm \
    pdbfixer \
    pdb2pqr \
    rdkit

# その他のパッケージ
RUN conda install -n in_silico -y \
    requests \
    urllib3 \
    tqdm \
    typing-extensions \
    jupyter \
    ipykernel \
    notebook \
    ipywidgets \
    plotly \
    bokeh \
    altair \
    sympy \
    networkx \
    lxml \
    h5py

# pipでcondaにないパッケージをインストール（Python 3.13対応版）
RUN /opt/conda/envs/in_silico/bin/pip install \
    py3Dmol==2.0.4

# download smina static version and add execution permission
RUN curl -Lo /usr/local/bin/smina https://sourceforge.net/projects/smina/files/smina.static/download
RUN chmod +x /usr/local/bin/smina

# Install mgltools
RUN conda install -n in_silico -c bioconda -y \
    mgltools

# Create a working directory
WORKDIR /workspace

CMD ["smina", "--version"]
# デフォルト動作（継続実行）
CMD ["tail", "-f", "/dev/null"]
