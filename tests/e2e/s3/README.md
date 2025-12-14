# S3 End-to-End Tests

This directory contains end-to-end tests for omnibenchmark's S3 remote storage functionality. The tests are designed to be minimal, fast, and focused on CLI signatures rather than comprehensive pipeline testing.

## Overview

The S3 e2e tests validate:
- CLI signature for `--use-remote-storage` flag
- Basic S3 operations (bucket creation, upload, download)
- Integration between omnibenchmark and S3-compatible storage

## Architecture

### Local Development (MinIO)
- Uses Docker Compose to run MinIO locally
- Replaces testcontainers for better control and reliability
- Automatic setup of test users and credentials
- Runs on `http://localhost:9000`

### CI/Production (Remote S3)
- Uses real S3 or S3-compatible services
- Timestamp-based bucket names to avoid conflicts
- Configured via environment variables
- Supports AWS S3, MinIO, and other S3-compatible services

## Quick Start

### Local Testing with MinIO

The easiest way to run S3 tests locally:

```bash
# From the omnibenchmark root directory
pixi run test-e2e-s3-local
```

Or manually:

```bash
cd tests/e2e/s3
./run-local-s3-test.sh
```

This will:
1. Check if MinIO is already running on port 9000
2. If not found, start MinIO containers with Docker Compose
3. Set up test users and credentials automatically
4. Run the S3 e2e tests
5. Clean up only containers it started (preserves existing MinIO)

#### Using Existing MinIO Container

If you already have MinIO running (e.g., from `../minio-data/run-minio.sh`):

```bash
# Set credentials if needed
export OB_STORAGE_S3_ACCESS_KEY=your-access-key
export OB_STORAGE_S3_SECRET_KEY=your-secret-key

# Run tests against existing MinIO
pixi run test-e2e-s3
```

The test script will detect your running MinIO and use it without starting new containers.

### Manual Setup

If you prefer to manage MinIO manually:

```bash
# Start MinIO
cd tests/e2e/s3
docker-compose up -d

# Wait for setup to complete
sleep 15

# Source credentials
source /tmp/minio-credentials

# Run tests
cd ../../../..  # back to project root
pixi run test-e2e-s3
```

## Remote S3 Testing

For testing against real S3 or other remote services:

1. Copy the example environment file:
   ```bash
   cp tests/e2e/s3/.env.example tests/e2e/s3/.env
   ```

2. Edit `.env` with your credentials:
   ```bash
   OB_E2E_USE_REMOTE_S3=true
   OB_STORAGE_S3_ACCESS_KEY=your-access-key
   OB_STORAGE_S3_SECRET_KEY=your-secret-key
   OB_STORAGE_S3_ENDPOINT_URL=https://s3.amazonaws.com
   ```

3. Source the environment and run tests:
   ```bash
   source tests/e2e/s3/.env
   pixi run test-e2e-s3
   ```

## Files

### Test Configuration
- `06_s3_minimal.yaml` - Minimal benchmark configuration for fast testing
- `test_06_s3_minimal.py` - Main test file with S3 e2e tests

### Docker Setup (Local Development)
- `docker-compose.yml` - MinIO container configuration
- `setup-minio.sh` - MinIO initialization script
- `run-local-s3-test.sh` - Convenience script for local testing

### Configuration
- `.env.example` - Example environment configuration
- `README.md` - This file

## Test Cases

All tests are marked with `@pytest.mark.e2e_s3` for separate execution from other e2e tests.

### `test_s3_remote_storage_cli_signature`
**Focus**: CLI signature validation
- Tests that `--use-remote-storage` flag works with S3 configuration
- Validates basic command execution and output creation
- Minimal pipeline to ensure speed

### `test_s3_remote_storage_without_flag_fails`
**Focus**: CLI validation
- Tests behavior when S3 config is used without `--use-remote-storage`
- Ensures proper error handling or graceful fallback

### `test_s3_bucket_creation_and_access`
**Focus**: S3 integration
- Tests bucket creation, file upload, and download operations
- Validates core S3 functionality beyond CLI signatures

## Environment Variables

### Test Control
- `OB_E2E_USE_REMOTE_S3` - Set to `true` to use remote S3 instead of local MinIO
- `OB_E2E_KEEP_FILES` - Set to `true` to keep test files for debugging
- `OB_E2E_BUCKET_PREFIX` - Prefix for CI bucket names (default: `obdata-ci`)
- `CI` - Set to `true` to automatically use remote S3 (for CI environments)

### S3 Configuration
- `OB_STORAGE_S3_ACCESS_KEY` - S3 access key
- `OB_STORAGE_S3_SECRET_KEY` - S3 secret key  
- `OB_STORAGE_S3_ENDPOINT_URL` - S3 endpoint URL (defaults to `eu-central-1` region)
- `AWS_DEFAULT_REGION` - AWS region (defaults to `eu-central-1`)

### GitHub Actions Context
- `GITHUB_PR_NUMBER` - PR number for bucket naming in CI
- `CI_MERGE_REQUEST_IID` - GitLab MR number (alternative to GitHub PR)

### Timeouts
- `OB_E2E_MINIO_STARTUP_TIMEOUT` - MinIO startup timeout in seconds (default: 120)
- `OB_E2E_TEST_TIMEOUT` - Test timeout in seconds (default: 300)

## CI Integration

The project includes a complete GitHub Actions workflow at `.github/workflows/s3-tests.yml` that:

1. **Runs on every PR** that touches S3-related code
2. **Validates S3 configuration** before running actual tests
3. **Runs S3 e2e tests** if credentials are configured
4. **Automatically cleans up** test buckets after completion
5. **Provides clear summaries** with setup instructions

#### Required Repository Secrets

Configure these secrets in your GitHub repository settings:

- `OB_STORAGE_S3_ACCESS_KEY` - S3 access key
- `OB_STORAGE_S3_SECRET_KEY` - S3 secret key  
- `OB_STORAGE_S3_ENDPOINT_URL` - S3 endpoint (optional, defaults to EU region)

#### Manual Workflow Trigger

You can also trigger the workflow manually with custom settings:

```bash
# Via GitHub UI: Actions → S3 E2E Tests → Run workflow
# Specify custom bucket prefix and debug options
```

### Bucket Management in CI

The tests automatically create unique bucket names:
- Format: `obdata-ci-prXX-YYYYMMDD-HHMMSS-sss` 
- Example: `obdata-ci-pr123-20241214-143022-456`

Features:
- **PR-specific**: Each PR gets its own bucket namespace
- **Timestamped**: No conflicts between test runs
- **Auto-cleanup**: Old buckets are automatically deleted after 1 hour

## Design Principles

### Minimal and Fast
- Uses smallest possible dataset (value: 50 vs 100+ in other tests)
- Single data module and single method module
- Focus on CLI signatures rather than complex pipelines

### Testcontainer Replacement
- Docker Compose provides better control than testcontainers
- Easier to debug and maintain
- More reliable container lifecycle management
- Explicit credential handling
- **Existing container support**: Works with already-running MinIO instances

### Environment Flexibility
- Same test code works for local MinIO and remote S3
- Configuration via environment variables
- Graceful fallbacks and clear error messages

### Functional Programming Approach
- Immutable configuration objects
- Pure functions for environment setup
- Clear separation of concerns
- Comprehensive error handling

## Troubleshooting

### MinIO Won't Start
```bash
# Check if port 9000 is already in use
lsof -i :9000

# Clean up any existing containers
docker-compose down -v --remove-orphans
```

### Credentials Not Working
```bash
# Check if credentials file was created
cat /tmp/minio-credentials

# Manually check MinIO is responding
curl http://localhost:9000/minio/health/live
```

### Tests Fail with Timeout
```bash
# Increase timeouts in .env
OB_E2E_MINIO_STARTUP_TIMEOUT=180
OB_E2E_TEST_TIMEOUT=600
```

### Debug Test Output
```bash
# Run with keep_files=True for debugging
OB_E2E_KEEP_FILES=true pixi run test-e2e-s3
```

## Running Tests

### Pixi Commands

The project provides several pixi tasks for different testing scenarios:

```bash
# Run only S3 e2e tests
pixi run test-e2e-s3

# Run S3 tests with local MinIO setup
pixi run test-e2e-s3-local

# Run all e2e tests EXCEPT S3 tests  
pixi run test-e2e

# Run ALL e2e tests including S3
pixi run test-e2e-all
```

### Pytest Markers

```bash
# Run only S3 e2e tests
pytest -m e2e_s3 -v

# Run all e2e tests except S3
pytest -m 'e2e and not e2e_s3' -v

# Run S3 config validation (no S3 setup required)
pytest tests/e2e/test_s3_config_validation.py -v
```

## Future Enhancements

1. **Testcontainer Removal**: Complete migration away from testcontainers ✅
2. **Multi-Backend Testing**: Support for different S3-compatible services
3. **Performance Benchmarking**: Add timing metrics for S3 operations
4. **Cleanup Jobs**: Automatic cleanup of old test buckets in CI ✅
5. **Integration Tests**: More comprehensive pipeline testing with S3