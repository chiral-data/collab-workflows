# Node 5: Run P2Rank pocket prediction
FROM eclipse-temurin:latest

# Install Python + required tooling
RUN apt-get update && \
    apt-get install -y --no-install-recommends python3 python3-pip python3-venv && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Create a virtual environment for Python packages
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir pandas

# Make the venv available in PATH
ENV PATH="/opt/venv/bin:$PATH"

WORKDIR /workspace
COPY . /workspace
RUN chmod +x /workspace/run.sh

# Create non-root user WITHOUT specifying UID
RUN useradd -m worker
USER worker

ENTRYPOINT ["bash", "/workspace/run.sh"]
