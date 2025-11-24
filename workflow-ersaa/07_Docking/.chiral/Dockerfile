# Node 7: Docking with AutoDock Vina
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip wget unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install pandas

# Install Vina
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    -O /usr/local/bin/vina && chmod +x /usr/local/bin/vina

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["bash", "run.sh"]
