FROM python:3.11-slim

LABEL maintainer="jlanej"
LABEL description="Cross-species contamination detection pipeline"

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
COPY csc/ csc/
RUN pip install --no-cache-dir . && \
    pip install --no-cache-dir pysam

# Install test dependencies (useful for CI; can be excluded in production)
COPY tests/ tests/
RUN pip install --no-cache-dir ".[test]"

# Verify installation
RUN csc-extract --version && samtools --version | head -1

# Verify all modules are importable
RUN python -c "import csc; import csc.extract; import csc.classify; import csc.aggregate; import csc.detect; import csc.utils; import csc.config"

ENTRYPOINT ["csc-extract"]
