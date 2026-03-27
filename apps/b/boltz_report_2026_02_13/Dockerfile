# Boltz/BoltzGen Report Generator
# Reused from container-images-for-potter/reports/boltz_report/Dockerfile_python
FROM python:3.8

# Install Python dependencies
COPY requirements.txt /workspace/requirements.txt
RUN pip install --no-cache-dir -r /workspace/requirements.txt

# Install system dependencies including jq for JSON parsing
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a working directory
WORKDIR /workspace

CMD ["tail", "-f", "/dev/null"]
