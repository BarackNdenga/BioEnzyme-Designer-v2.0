#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# BioEnzyme Designer v2.0 — Google Cloud Deployment Script
# Deploy to Google Cloud Run (serverless container)
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────
PROJECT="${GCP_PROJECT:-}"
REGION="${GCP_REGION:-us-central1}"
SERVICE_NAME="bioenzyme"
IMAGE="gcr.io/${PROJECT}/${SERVICE_NAME}:v2.0"
MEMORY="2Gi"
CPU="2"
CONCURRENCY=10
MAX_INSTANCES=5

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  BioEnzyme Designer v2.0 — Google Cloud Run Deployment        ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

# ── Step 1: Build and push Docker image ───────────────────────────────────
echo "[1/4] Building Docker image..."
docker build -t "${IMAGE}" \
  -f cloud/Dockerfile .

echo "[2/4] Pushing image to GCR..."
docker push "${IMAGE}"
echo "  → Image: ${IMAGE}"

# ── Step 2: Deploy to Cloud Run ───────────────────────────────────────────
echo "[3/4] Deploying to Cloud Run..."

gcloud run deploy "${SERVICE_NAME}" \
  --image "${IMAGE}" \
  --region "${REGION}" \
  --memory "${MEMORY}" \
  --cpu "${CPU}" \
  --concurrency "${CONCURRENCY}" \
  --max-instances "${MAX_INSTANCES}" \
  --timeout "600" \
  --port 8000 \
  --allow-unauthenticated \
  --add-cloudsql-instances "" \
  --set-env-vars="LOG_LEVEL=info,BRENDA_API_KEY=${BRENDA_API_KEY:-}" \
  --labels="app=bioenzyme,version=2.0"

# ── Step 3: Print service URL ─────────────────────────────────────────────
echo ""
echo "[4/4] Getting service URL..."
SERVICE_URL=$(gcloud run services describe "${SERVICE_NAME}" \
  --region "${REGION}" \
  --format="value(status.url)")

echo ""
echo "✅ Deployment complete!"
echo "   Service URL: ${SERVICE_URL}"
echo "   Health check: ${SERVICE_URL}/health"
echo "   API docs:     ${SERVICE_URL}/docs"
echo ""
echo "Example usage:"
echo "  curl -X POST ${SERVICE_URL}/analyze \\"
echo "    -F 'file=@enzyme.pdb' \\"
echo "    -d 'improve=activity' \\"
echo "    -d 'data_source=combined'"
