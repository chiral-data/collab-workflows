# Use the official TensorFlow GPU image
FROM tensorflow/tensorflow:2.20.0-gpu

# Set the working directory in the container
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir numpy pandas rdkit scikit-learn


# Create a working directory
WORKDIR /workspace

CMD ["python", "-c", "import tensorflow as tf; print(tf.__version__)"]
