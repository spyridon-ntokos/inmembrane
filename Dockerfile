# Minimal Dockerfile for the modern inmembrane fork (SerraPHIM integration)
# -------------------------------------------------------------------------
# This image:
#   - Provides Python 3 + HMMER
#   - Installs inmembrane and its Python dependencies
#
# It does NOT bundle SignalP 6.0, TMbed or DeepLocPro, because those have
# their own installation and/or license requirements. You are expected to:
#   - Either mount a host environment where these tools are available
#   - Or build your own derived image that installs them as needed.

FROM python:3.10-slim

# System deps
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        hmmer \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# App directory
WORKDIR /opt/inmembrane
COPY . .

# Python deps + inmembrane itself
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt && \
    pip install .

# Default working dir where users can mount data
WORKDIR /data

# CLI entrypoint
ENTRYPOINT ["inmembrane_scan"]
