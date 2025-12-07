
# Use an official Python runtime as a parent image
FROM python:3.11-slim-buster

# Set the working directory in the container
WORKDIR /app

# Install pocketeer and its dependencies
RUN pip install pocketeer[vis]

# Copy the current directory contents into the container at /app
COPY complete_notebook_example.py .

# Command to run the Python script
CMD ["python", "complete_notebook_example.py"]
