# ═══════════════════════════════════════════════════════════════════════════
# BioEnzyme Designer v2.0 — Dockerfile
# Multi-stage build for cloud deployment (AWS ECS, Google Cloud Run, etc.)
# ═══════════════════════════════════════════════════════════════════════════

FROM python:3.11-slim AS base

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libxml2-dev \
    libxslt-dev \
    zlib1g-dev \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Install the package in editable mode
RUN pip install -e .

# Expose port for API
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Default command: run the API server
CMD ["python", "-m", "bioenzyme_v2.api.rest_api", "--host", "0.0.0.0", "--port", "8000"]


# ═══════════════════════════════════════════════════════════════════════════
# Alternative: lightweight version for CLI-only use
# ═══════════════════════════════════════════════════════════════════════════
FROM python:3.11-slim AS cli

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements-cli.txt .
RUN pip install --no-cache-dir -r requirements-cli.txt

COPY . .
RUN pip install -e .

ENTRYPOINT ["python", "-m", "bioenzyme_v2.cli"]
