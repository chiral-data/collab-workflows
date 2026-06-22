FROM ubuntu:22.04

# Java runtime (required by P2Rank) + Python 3 (used by run.sh to select best model)
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre-headless \
    wget \
    python3 \
    && rm -rf /var/lib/apt/lists/*

# P2Rank 2.4.2 — baked in at build time, not downloaded at runtime
RUN wget -q https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz \
        -O /tmp/p2rank.tar.gz \
    && tar -xzf /tmp/p2rank.tar.gz -C /opt \
    && rm /tmp/p2rank.tar.gz \
    && chmod +x /opt/p2rank_2.4.2/prank

ENV PATH="/opt/p2rank_2.4.2:$PATH"

WORKDIR /workspace
