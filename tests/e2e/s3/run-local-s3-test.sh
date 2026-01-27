#!/usr/bin/env bash
set -e

# Convenience script for running S3 e2e tests locally with MinIO
# This script handles the setup and teardown of the local MinIO environment

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

echo "=== omnibenchmark S3 E2E Test Runner ==="
echo "Project root: $PROJECT_ROOT"
echo "S3 test dir: $SCRIPT_DIR"

# Function to cleanup on exit
cleanup() {
    echo ""
    echo "=== Cleaning up ==="
    
    # Only cleanup containers we started (if STARTED_CONTAINERS is set)
    if [ "$STARTED_CONTAINERS" = "true" ]; then
        echo "Stopping containers we started..."
        cd "$SCRIPT_DIR"
        docker compose down -v --remove-orphans || true
    else
        echo "Leaving existing MinIO containers running..."
    fi
    
    # Clean up any leftover environment variables
    unset OB_STORAGE_S3_ACCESS_KEY
    unset OB_STORAGE_S3_SECRET_KEY
    unset OB_STORAGE_S3_ENDPOINT_URL
    unset STARTED_CONTAINERS
    
    echo "Cleanup complete."
}

# Set trap for cleanup
trap cleanup EXIT INT TERM

# Check prerequisites
echo ""
echo "=== Checking Prerequisites ==="

if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is required but not installed."
    exit 1
fi

if ! command -v docker &> /dev/null || ! docker compose version &> /dev/null; then
    echo "WARNING: Docker Compose not found - will only work with existing MinIO"
fi

if ! command -v nc &> /dev/null; then
    echo "WARNING: netcat (nc) not found - container detection may be unreliable"
fi

echo "✓ Docker and Docker Compose are available"

# Check if MinIO is already running
echo ""
echo "=== Checking for MinIO ==="

# Check if MinIO is running on port 9000
if nc -z localhost 9000 2>/dev/null; then
    echo "✓ Found existing MinIO on port 9000"
    
    # Try to use existing credentials from environment or file
    if [ -n "$OB_STORAGE_S3_ACCESS_KEY" ] && [ -n "$OB_STORAGE_S3_SECRET_KEY" ]; then
        echo "✓ Using credentials from environment variables"
        export OB_STORAGE_S3_ENDPOINT_URL="http://localhost:9000"
    elif [ -f "/tmp/minio-credentials" ]; then
        echo "✓ Using credentials from file"
        source /tmp/minio-credentials
        export OB_STORAGE_S3_ACCESS_KEY
        export OB_STORAGE_S3_SECRET_KEY
        export OB_STORAGE_S3_ENDPOINT_URL
    else
        echo "⚠ No credentials found, using MinIO defaults"
        export OB_STORAGE_S3_ACCESS_KEY="minioadmin"
        export OB_STORAGE_S3_SECRET_KEY="minioadmin123"
        export OB_STORAGE_S3_ENDPOINT_URL="http://localhost:9000"
    fi
    
    echo "  Access Key: $OB_STORAGE_S3_ACCESS_KEY"
    echo "  Endpoint: $OB_STORAGE_S3_ENDPOINT_URL"
    
else
    echo "No MinIO found on port 9000, starting new containers..."
    cd "$SCRIPT_DIR"

    # Stop any existing containers
    docker compose down -v --remove-orphans || true

    # Start fresh containers
    echo "Starting MinIO containers..."
    docker compose up -d
    STARTED_CONTAINERS="true"

    # Wait for setup to complete
    echo "Waiting for MinIO setup to complete..."
    sleep 15

    # Check if credentials file was created
    if [ -f "/tmp/minio-credentials" ]; then
        echo "✓ MinIO credentials created successfully"
        source /tmp/minio-credentials
        export OB_STORAGE_S3_ACCESS_KEY
        export OB_STORAGE_S3_SECRET_KEY
        export OB_STORAGE_S3_ENDPOINT_URL
        echo "  Access Key: $OB_STORAGE_S3_ACCESS_KEY"
        echo "  Endpoint: $OB_STORAGE_S3_ENDPOINT_URL"
    else
        echo "⚠ Credentials file not found, using defaults"
        export OB_STORAGE_S3_ACCESS_KEY="minioadmin"
        export OB_STORAGE_S3_SECRET_KEY="minioadmin123"
        export OB_STORAGE_S3_ENDPOINT_URL="http://localhost:9000"
    fi
fi

# Run the tests
echo ""
echo "=== Running S3 E2E Tests ==="
cd "$PROJECT_ROOT"

# Set environment to use local MinIO
export OB_E2E_USE_REMOTE_S3=false

# Check if pixi is available and use it, otherwise use pytest directly
if command -v pixi &> /dev/null; then
    echo "Using pixi to run tests..."
    pixi run pytest -m e2e_s3 -v -s
else
    echo "Using pytest directly..."
    pytest -m e2e_s3 -v -s
fi

echo ""
echo "=== S3 E2E Tests Complete ==="
echo "MinIO Console: http://localhost:9001"
echo "Username: minioadmin"  
echo "Password: minioadmin123"
echo ""
if [ "$STARTED_CONTAINERS" = "true" ]; then
    echo "The containers will be cleaned up automatically when this script exits."
else
    echo "Using existing MinIO containers - they will remain running."
fi