FROM chiral/vina-base:1.0
SHELL ["/bin/bash","-c"]

# Switch to working directory BEFORE copy
WORKDIR /workspace

# Copy all files while still root
COPY --chmod=755 run.sh \workspace\run.sh

USER worker

# Final entrypoint
ENTRYPOINT ["bash", "/workspace/run.sh"]
