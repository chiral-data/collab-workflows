# Node 8: Rank Vina results
FROM python:3.10-slim

RUN pip install pandas openpyxl

RUN useradd -m -u 1000 worker
USER worker

WORKDIR /workspace
COPY . /workspace

RUN chmod +x run.sh

ENTRYPOINT ["bash", "run.sh"]
