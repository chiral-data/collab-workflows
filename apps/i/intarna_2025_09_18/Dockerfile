# Use the official IntaRNA biocontainer from Quay.io
# This provides a fully tested and optimized IntaRNA installation
FROM quay.io/biocontainers/intarna:3.4.1--pl5321h077b44d_3

# Set working directory
WORKDIR /workspace

# The biocontainer already includes IntaRNA and all dependencies
# IntaRNA is available in the PATH

# Python is available for JSON processing (no additional packages needed)

# Default command, used for verification
CMD ["IntaRNA", "--help"]