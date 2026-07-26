#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# BioEnzyme Designer v2.0 — AWS Deployment Script
# Deploy to AWS ECS Fargate behind an Application Load Balancer
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────
PROJECT="bioenzyme"
REGION="${AWS_REGION:-us-east-1}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PROJECT}"
IMAGE_TAG="v2.0"
ECS_CLUSTER="${PROJECT}-cluster"
ECS_SERVICE="${PROJECT}-service"
TASK_DEF="${PROJECT}-task"
SUBNETS="${AWS_SUBNETS:-}"
SECURITY_GROUP="${AWS_SG:-}"
LOAD_BALANCER="${AWS_LB_ARN:-}"
LISTENER_PORT=8000

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  BioEnzyme Designer v2.0 — AWS ECS Fargate Deployment         ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

# ── Step 1: Build and push Docker image ───────────────────────────────────
echo "[1/6] Building Docker image..."
docker build -t "${ECR_REPO}:${IMAGE_TAG}" \
  --build-arg TARGET=base \
  -f cloud/Dockerfile .

echo "[2/6] Logging into ECR..."
aws ecr get-login-password --region "${REGION}" | \
  docker login --username AWS --password-stdin "${ECR_REPO}"

echo "[3/6] Pushing image to ECR..."
docker push "${ECR_REPO}:${IMAGE_TAG}"
echo "  → Image: ${ECR_REPO}:${IMAGE_TAG}"

# ── Step 2: Create ECS cluster (if not exists) ────────────────────────────
echo "[4/6] Setting up ECS cluster..."
if ! aws ecs describe-clusters --clusters "${ECS_CLUSTER}" --region "${REGION}" \
     >/dev/null 2>&1; then
    aws ecs create-cluster \
      --cluster-name "${ECS_CLUSTER}" \
      --capacity-providers FARGATE FARGATE_SPOT \
      --region "${REGION}"
    echo "  → Cluster created: ${ECS_CLUSTER}"
else
    echo "  → Cluster exists: ${ECS_CLUSTER}"
fi

# ── Step 3: Register task definition ──────────────────────────────────────
echo "[5/6] Registering task definition..."
cat > /tmp/task-def.json <<EOF
{
  "family": "${TASK_DEF}",
  "networkMode": "awsvpc",
  "requiresCompatibilities": ["FARGATE"],
  "cpu": "1024",
  "memory": "2048",
  "executionRoleArn": "arn:aws:iam::${ACCOUNT_ID}:role/ecsTaskExecutionRole",
  "containerDefinitions": [
    {
      "name": "bioenzyme-api",
      "image": "${ECR_REPO}:${IMAGE_TAG}",
      "portMappings": [
        {
          "containerPort": 8000,
          "protocol": "tcp"
        }
      ],
      "essential": true,
      "healthCheck": {
        "command": ["CMD-SHELL", "curl -f http://localhost:8000/health || exit 1"],
        "interval": 30,
        "timeout": 10,
        "retries": 3
      },
      "environment": [
        {"name": "LOG_LEVEL", "value": "info"}
      ],
      "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-group": "/ecs/${PROJECT}",
          "awslogs-region": "${REGION}",
          "awslogs-stream-prefix": "ecs"
        }
      }
    }
  ]
}
EOF

TASK_DEF_ARN=$(aws ecs register-task-definition \
  --cli-input-json file:///tmp/task-def.json \
  --region "${REGION}" \
  --query "taskDefinition.taskDefinitionArn" \
  --output text)
echo "  → Task definition: ${TASK_DEF_ARN}"

# ── Step 4: Create/update ECS service ─────────────────────────────────────
echo "[6/6] Deploying ECS service..."
if aws ecs describe-services \
     --cluster "${ECS_CLUSTER}" \
     --services "${ECS_SERVICE}" \
     --region "${REGION}" >/dev/null 2>&1; then
    aws ecs update-service \
      --cluster "${ECS_CLUSTER}" \
      --service "${ECS_SERVICE}" \
      --task-definition "${TASK_DEF_ARN}" \
      --desired-count 1 \
      --region "${REGION}"
    echo "  → Service updated: ${ECS_SERVICE}"
else
    aws ecs create-service \
      --cluster "${ECS_CLUSTER}" \
      --service-name "${ECS_SERVICE}" \
      --task-definition "${TASK_DEF_ARN}" \
      --desired-count 1 \
      --launch-type FARGATE \
      --network-configuration "awsvpcConfiguration={subnets=[${SUBNETS}],securityGroups=[${SECURITY_GROUP}],assignPublicIp=ENABLED}" \
      --load-balancers "targetGroupArn=${LOAD_BALANCER},containerName=bioenzyme-api,containerPort=8000" \
      --region "${REGION}"
    echo "  → Service created: ${ECS_SERVICE}"
fi

echo ""
echo "✅ Deployment complete!"
echo "   Wait for the service to stabilise (~2 min)."
echo "   Check: aws ecs describe-services --cluster ${ECS_CLUSTER} --services ${ECS_SERVICE} --region ${REGION}"
