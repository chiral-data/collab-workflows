# Use the official CopraRNA biocontainer from Quay.io
# This provides a fully tested CopraRNA installation with all dependencies
# Version 2.1.3 is the latest stable release
FROM quay.io/biocontainers/coprarna:2.1.3--6

# Set working directory
WORKDIR /workspace

# The biocontainer already includes:
# - CopraRNA 2.1.3
# - IntaRNA (>2.2)
# - BLAST, ClustalO, MAFFT
# - Perl with required modules
# - R with required packages
# - All other dependencies

# Default command to show help
CMD ["CopraRNA2.pl"]