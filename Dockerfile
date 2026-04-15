FROM python:3.12-slim

LABEL maintainer="jlanej"
LABEL description="Cross-species contamination detection pipeline"

ARG KRAKEN2_VERSION=v2.1.3

# Install samtools, Kraken2 build dependencies, and runtime libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        perl \
        zlib1g-dev \
        samtools \
    && git clone --depth 1 --branch "${KRAKEN2_VERSION}" \
        https://github.com/DerrickWood/kraken2.git /tmp/kraken2 \
    && cd /tmp/kraken2 && ./install_kraken2.sh /usr/local/bin \
    && rm -rf /tmp/kraken2 \
    && rm -rf /var/lib/apt/lists/*

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
RUN csc-extract --version && csc-classify --version && \
    samtools --version | head -1 && \
    kraken2 --version | head -1

# Verify all modules are importable
RUN python -c "import csc; import csc.extract; import csc.classify; import csc.aggregate; import csc.detect; import csc.utils; import csc.config"

# Kraken2 databases should be mounted as a volume at runtime:
#   docker run -v /path/to/kraken2_db:/data/kraken2_db ...
VOLUME ["/data/kraken2_db"]

ENTRYPOINT []
CMD ["/bin/bash"]
