#!/bin/sh
set -e

echo "Setting up RustFS credentials for omnibenchmark e2e tests..."

# Save credentials to a file that can be read by tests
cat > /tmp/rustfs-credentials <<EOF
OB_STORAGE_S3_ACCESS_KEY=rustfsadmin
OB_STORAGE_S3_SECRET_KEY=rustfsadmin
OB_STORAGE_S3_ENDPOINT_URL=http://localhost:9000
EOF

echo "RustFS setup complete!"
echo "Credentials saved to /tmp/rustfs-credentials"
echo "Access Key: rustfsadmin"
echo "Secret Key: rustfsadmin"
echo "Endpoint: http://localhost:9000"
echo "Console: http://localhost:9001"
