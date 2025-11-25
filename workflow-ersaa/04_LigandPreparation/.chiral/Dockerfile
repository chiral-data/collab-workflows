FROM chiral/vina-base:1.0
SHELL ["/bin/bash","-c"]

USER worker
WORKDIR /workspace

COPY . /workspace
RUN chmod +x /workspace/run.sh

ENTRYPOINT ["bash","/workspace/run.sh"]
