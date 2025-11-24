# Node 5: Run P2Rank pocket prediction
FROM openjdk:17-slim

RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install pandas

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["bash", "run.sh"]
