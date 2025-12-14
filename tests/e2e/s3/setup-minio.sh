#!/bin/sh
set -e

echo "Setting up MinIO for omnibenchmark e2e tests..."

# Wait for MinIO to be ready
until mc ready local; do
  echo "Waiting for MinIO to be ready..."
  sleep 2
done

echo "MinIO is ready!"

# Enable versioning support globally for better S3 compatibility
echo "Configuring MinIO settings..."
mc admin config set local api requests_max=10000
mc admin service restart local

# Wait for restart
sleep 5
until mc ready local; do
  echo "Waiting for MinIO restart..."
  sleep 2
done

# Use root credentials directly for tests (simplifies permissions and cleanup)
echo "Using root credentials for tests..."

# Save credentials to a file that can be read by tests
cat > /tmp/minio-credentials <<EOF
OB_STORAGE_S3_ACCESS_KEY=minioadmin
OB_STORAGE_S3_SECRET_KEY=minioadmin123
OB_STORAGE_S3_ENDPOINT_URL=http://localhost:9000
EOF

echo "MinIO setup complete!"
echo "Credentials saved to /tmp/minio-credentials"
echo "Access Key: minioadmin (root)"
echo "Secret Key: minioadmin123 (root)" 
echo "Endpoint: http://localhost:9000"
echo "Console: http://localhost:9001"
echo "Note: Using root credentials for full S3 API compatibility"