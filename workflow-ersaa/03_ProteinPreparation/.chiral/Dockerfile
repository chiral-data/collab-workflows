# rely on the heavy base image that already contains MGLTools/OpenMM/OpenBabel
FROM chiral/vina-base:1.0
SHELL ["/bin/bash","-c"]

# create a light activation wrapper and keep non-root
USER worker
WORKDIR /workspace

COPY . /workspace
RUN chmod +x /workspace/run.sh

ENTRYPOINT ["bash","/workspace/run.sh"]
