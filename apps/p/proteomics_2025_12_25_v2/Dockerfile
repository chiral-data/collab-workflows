FROM python:3.11-slim

RUN pip install --no-cache-dir \
    pandas \
    numpy \
    scikit-learn \
    scipy \
    seaborn \
    matplotlib \
    scikit-posthocs

WORKDIR /workspace

CMD ["python", "-c", "import pandas; import sklearn; import scipy; print('proteomics ready')"]
