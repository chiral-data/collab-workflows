FROM python:3.11-slim-bookworm

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies (excluding streamlit-related packages)
RUN pip install --no-cache-dir \
    plotly \
    pandas \
    numpy \
    scipy \
    Pillow \
    matplotlib \
    reportlab \
    requests \
    opencv-python-headless \
    cellpose \
    torch \
    torchvision \
    scikit-image \
    einops

# Create a working directory
WORKDIR /workspace

CMD ["tail", "-f", "/dev/null"]

