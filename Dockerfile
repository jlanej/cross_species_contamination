FROM python:3.11-slim

LABEL maintainer="jlanej"
LABEL description="Unmapped-read extraction pipeline for cross-species contamination detection"

# Install samtools and dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools \
        libhts-dev \
        gcc \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python package
COPY pyproject.toml README.md LICENSE ./
COPY extract_unmapped/ extract_unmapped/
RUN pip install --no-cache-dir . && \
    pip install --no-cache-dir pysam

# Install test dependencies (useful for CI; can be excluded in production)
COPY tests/ tests/
RUN pip install --no-cache-dir ".[test]"

# Verify installation
RUN extract-unmapped --version && samtools --version | head -1

ENTRYPOINT ["extract-unmapped"]
