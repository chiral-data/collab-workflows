# Use a lightweight Python base image
FROM python:3.12-slim-bookworm

# Set the working directory in the container
WORKDIR /app

# Install system dependencies that RDKit might need
# For example, `libsm6` and `libxrender1` are often required for RDKit with matplotlib (though not directly used here, it's good practice for image processing libraries)
RUN apt-get update && apt-get install -y \
    build-essential \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir numpy pandas rdkit tensorflow scikit-learn


# Create a working directory
WORKDIR /workspace

CMD ["python", "-c", "import tensorflow as tf; print(tf.__version__)"]
