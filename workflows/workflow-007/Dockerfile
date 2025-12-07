# Use an official Python runtime as a parent image
FROM python:3.11-slim-bookworm

# Set the working directory in the container
WORKDIR /app

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

# Copy the current directory contents into the container at /app
COPY complete_notebook_example.py .

# Command to run the Python script
CMD ["python", "complete_notebook_example.py"]
