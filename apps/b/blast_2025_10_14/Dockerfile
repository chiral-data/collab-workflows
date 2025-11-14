# Use the official NCBI BLAST image
# This is the authoritative BLAST image maintained by NCBI
# Version 2.17.0 from July 2025
FROM ncbi/blast:2.17.0

# Install python3 and jq for JSON config parsing
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# The image already includes the complete BLAST+ suite:
# blastn, blastp, blastx, tblastn, tblastx, makeblastdb, blastdbcmd, etc.

# Default command, used for verification
CMD ["blastn", "-version"]
