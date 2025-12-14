# S3 E2E Tests - Quick Usage Guide

This guide shows you how to run the S3 end-to-end tests in different scenarios.

## 🚀 Quick Start

### Option 1: Auto Setup with Local MinIO (Recommended)
```bash
pixi run test-e2e-s3-local
```
This handles everything automatically - starts MinIO if needed, runs tests, cleans up.

### Option 2: Use Existing MinIO Container
If you already have MinIO running (e.g., from `../minio-data/run-minio.sh`):
```bash
pixi run test-e2e-s3
```

### Option 3: Use Remote S3
```bash
export OB_E2E_USE_REMOTE_S3=true
export OB_STORAGE_S3_ACCESS_KEY=your-key
export OB_STORAGE_S3_SECRET_KEY=your-secret
export OB_STORAGE_S3_ENDPOINT_URL=https://s3.eu-central-1.amazonaws.com
pixi run test-e2e-s3
```

## 🎯 What These Tests Do

The S3 e2e tests are **minimal and fast** - they focus on CLI signatures rather than comprehensive pipeline testing:

1. **CLI Signature Test**: Verifies `ob run benchmark -b config.yaml --use-remote-storage` works
2. **Error Handling Test**: Checks behavior when S3 config is used without `--use-remote-storage`
3. **S3 Integration Test**: Validates bucket creation, upload, and download operations

**Expected runtime**: < 2 minutes
**Expected outputs**: 2+ JSON files demonstrating S3 upload/download worked

## 📋 Prerequisites

### For Local MinIO
- Docker and Docker Compose
- Port 9000 available (or existing MinIO on that port)

### For Remote S3
- S3 credentials (access key + secret key)
- S3 endpoint URL (defaults to EU region)

## 🛠️ Advanced Usage

### Run with Different Test Configurations

```bash
# Keep test files for debugging
OB_E2E_KEEP_FILES=true pixi run test-e2e-s3

# Use custom bucket prefix
OB_E2E_BUCKET_PREFIX=my-test pixi run test-e2e-s3

# Specify different MinIO startup timeout
OB_E2E_MINIO_STARTUP_TIMEOUT=180 pixi run test-e2e-s3-local
```

### Run Individual Test Functions

```bash
# Test only CLI signature
pytest tests/e2e/test_06_s3_minimal.py::test_s3_remote_storage_cli_signature -v

# Test only error handling
pytest tests/e2e/test_06_s3_minimal.py::test_s3_remote_storage_without_flag_fails -v

# Test only S3 integration
pytest tests/e2e/test_06_s3_minimal.py::test_s3_bucket_creation_and_access -v
```

### Run with Different Pytest Markers

```bash
# Run all S3 e2e tests
pytest -m e2e_s3 -v

# Run all e2e tests EXCEPT S3
pytest -m 'e2e and not e2e_s3' -v

# Run config validation (no S3 setup required)
pytest tests/e2e/test_s3_config_validation.py -v
```

## 🔧 Troubleshooting

### MinIO Won't Start
```bash
# Check port usage
lsof -i :9000

# Clean up containers
cd tests/e2e/s3
docker-compose down -v --remove-orphans
```

### Credentials Issues
```bash
# Check MinIO credentials file
cat /tmp/minio-credentials

# Test MinIO connection manually
curl http://localhost:9000/minio/health/live
```

### Test Failures
```bash
# Run with debug output
OB_E2E_KEEP_FILES=true pixi run test-e2e-s3 -s -v

# Check container logs
docker-compose logs minio
```

## 🏗️ CI/CD Integration

### GitHub Actions
The project includes `.github/workflows/s3-tests.yml` that:
- Runs on every PR touching S3-related code
- Uses repository secrets for S3 credentials
- Creates unique buckets per PR: `obdata-ci-pr123-timestamp`
- Auto-cleans up test buckets

### Required Repository Secrets
- `OB_STORAGE_S3_ACCESS_KEY`
- `OB_STORAGE_S3_SECRET_KEY`
- `OB_STORAGE_S3_ENDPOINT_URL` (optional)

## 📁 File Structure

```
tests/e2e/s3/
├── docker-compose.yml     # MinIO container setup
├── setup-minio.sh         # MinIO initialization script
├── run-local-s3-test.sh   # Convenience script for local testing
├── .env.example           # Environment configuration template
├── README.md              # Detailed documentation
└── USAGE.md               # This file

tests/e2e/
├── test_06_s3_minimal.py       # Main S3 e2e tests
├── test_s3_config_validation.py # Config validation (no S3 required)
└── configs/
    ├── 06_s3_minimal.yaml          # Minimal S3 test configuration
    └── 06_s3_minimal.expected.json # Expected test results
```

## 🔍 What Gets Tested

### Command Line Interface
```bash
# This exact command is tested:
ob run benchmark -b 06_s3_minimal.yaml --use-remote-storage
```

### S3 Operations
- Bucket creation with unique names
- File upload during data stage
- File download during method stage
- Proper error handling

### Configuration Validation
- YAML structure and required fields
- S3 endpoint and bucket name format
- Software environment setup
- Expected output file patterns

## 💡 Tips

1. **Fast Iteration**: Use `pixi run test-e2e-s3` against existing MinIO for fastest development
2. **Debug Mode**: Always use `OB_E2E_KEEP_FILES=true` when investigating failures
3. **CI Testing**: The workflow automatically handles bucket cleanup and provides clear summaries
4. **Existing Containers**: The scripts detect and work with existing MinIO instances
5. **European Region**: Default S3 endpoint uses `eu-central-1` region

## 📚 Further Reading

- `README.md` - Complete documentation with architecture details
- `../../../06_s3_remote.py` - Original S3 configuration this is based on
- `.github/workflows/s3-tests.yml` - CI workflow configuration