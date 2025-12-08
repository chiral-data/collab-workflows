# Use an official Python runtime as a parent image
FROM python:3.11-slim-bookworm

# Install system dependencies for PyMOL
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1-mesa-glx \
    libglew2.2 \
    libpng16-16 \
    libfreetype6 \
    libxml2 \
    libglu1-mesa \
    && rm -rf /var/lib/apt/lists/*

# Install pocketeer, PyMOL, and imageio
RUN pip install --no-cache-dir pocketeer[vis] pymol-open-source imageio

# Create a working directory
WORKDIR /workspace

# Command to show pocketeer version
CMD ["python", "-c", "import pocketeer; print(pocketeer.__version__)"]
