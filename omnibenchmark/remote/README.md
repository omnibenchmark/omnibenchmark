# Remote Storage Module

The `omnibenchmark.remote` module provides functionality for managing remote storage (MinIO/S3), versioning benchmark results, archiving benchmarks, and handling file operations.

## Table of Contents

- [Overview](#overview)
- [Module Responsibilities](#module-responsibilities)
- [Remote Storage Operations](#remote-storage-operations)
  - [Version Management](#version-management)
  - [File Operations](#file-operations)
  - [Policy Generation](#policy-generation)
- [Archive Operations](#archive-operations)
- [Usage Examples](#usage-examples)
- [Configuration](#configuration)
- [Creating Access Keys](#creating-access-keys)

## Overview

This module handles two main concerns:

1. **Remote Storage Management**: Version-controlled S3/MinIO storage for benchmark results
2. **Archive Creation**: Packaging benchmarks with code, software dependencies, and results

## Module Responsibilities

### Remote Storage (`omnibenchmark.remote`)

- **Storage Abstraction**: Unified interface for S3/MinIO backends
- **Benchmark Versioning**: Create immutable snapshots with object tagging and retention policies
- **File Management**: Upload, download, list files with integrity checking (MD5 checksums)
- **Policy Management**: Generate least-privilege IAM policies for benchmark access
- **Credential Handling**: Secure access key management from environment variables or config files

### Archive (`omnibenchmark.remote.archive`)

- **Code Collection**: Clone Git repositories at specific commits
- **Software Collection**: Gather software environment definitions (conda, apptainer, envmodules)
- **Results Collection**: Download benchmark results from remote storage
- **Packaging**: Create compressed archives (ZIP, BZIP2, LZMA) of complete benchmarks

## Remote Storage Operations

### Version Management

Remote storage uses **version-controlled S3 buckets** where files are tagged with benchmark versions and protected with retention policies.

#### How Versioning Works

For each benchmark, a version-controlled S3 bucket is created. When you upload files, S3 assigns a `versionID` to each object:

```
/
├── out/
│   ├── f1.txt  (versionID: V2)
│   │           (versionID: V1)
│   └── f2.txt  (versionID: V2)
│               (versionID: V1)
```

Creating a benchmark version (snapshot) tags the newest file versions and write-protects them:

```
/
├── versions/
│   └── 0.1.csv              # Version manifest
└── out/
    ├── f1.txt
    │   (versionID: V2, tags: version=0.1, Retention: Governance)
    │   (versionID: V1)
    └── f2.txt
        (versionID: V2, tags: version=0.1, Retention: Governance)
        (versionID: V1)
```

The `0.1.csv` manifest contains: `name`, `version_id`, `last_modified`, `size`, `etag`.

#### Version Protection

Tagged versions are protected using **S3 Object Lock (Governance Mode)**:

- Tagged objects cannot be deleted or modified without special permissions
- Prevents accidental data loss
- Ensures reproducibility of benchmark results
- Allows rollback if version creation fails

### File Operations

The module provides three main file operations:

#### List Files

List available files for a benchmark:

```bash
ob remote files list -b benchmark.yaml
```

#### Download Files

Download files with checksum verification:

```bash
ob remote files download -b benchmark.yaml [--overwrite]
```

#### Checksum Validation

Verify file integrity using MD5 checksums:

```bash
ob remote files checksum -b benchmark.yaml
```

### Policy Generation

Generate **least-privilege IAM policies** for benchmark access.

**Generated Policy Capabilities:**
- ✅ Full S3 operations on benchmark bucket (read, write, list)
- ❌ Cannot bypass governance retention
- ❌ Cannot delete object versions or tags
- ❌ Cannot delete bucket or bucket policy

**Note**: Policy generation is an administrative operation. For production use, consider moving this to `obadmin` tooling.

## Archive Operations

Archives package complete benchmarks for distribution or long-term storage.

### Archive Components

- **Config**: Benchmark definition YAML file
- **Code**: Git repositories cloned at specific commits
- **Software**: 
  - EasyConfig files (for envmodules)
  - Conda environment YAML files
  - Apptainer container files (.sif)
- **Results**: Benchmark output files (from local or remote storage)

### Compression Options

- `none` - No compression (ZIP_STORED)
- `deflated` - Standard ZIP compression (ZIP_DEFLATED)
- `bzip2` - Better compression, slower (ZIP_BZIP2)
- `lzma` - Best compression, slowest (ZIP_LZMA)

## Usage Examples

### Version Management Commands

#### Create Benchmark Version

Create a new version snapshot of your benchmark results:

```bash
ob remote version create -b benchmark.yaml
```

This command:
- Tags the newest versions of all output files with the benchmark version
- Applies Governance retention policy to prevent deletion
- Creates a version manifest CSV file

#### List Available Versions

```bash
ob remote version list -b benchmark.yaml
```

Example output:
```
Available versions of benchmark.yaml:
     0.1
     0.2
     0.3
```

#### Compare Versions

Show differences between two benchmark versions:

```bash
ob remote version diff -b benchmark.yaml --version1 0.1 --version2 0.2
```

This displays a unified diff showing:
- Files added or removed between versions
- Changes in file sizes
- Modification timestamps

### File Management Commands

#### List Files

List all benchmark output files:

```bash
ob remote files list -b benchmark.yaml
```

Example output:
```
d41d8cd98f00b204e9800998ecf8427e out/data/D1/default/D1.txt.gz
a1b2c3d4e5f6g7h8i9j0k1l2m3n4o5p6 out/methods/M1/default/result.csv
```

#### Download Files

Download all benchmark files to local directory:

```bash
ob remote files download -b benchmark.yaml
```

Overwrite existing local files:

```bash
ob remote files download -b benchmark.yaml --overwrite
```

#### Verify Checksums

Generate and verify MD5 checksums for all benchmark outputs:

```bash
ob remote files checksum -b benchmark.yaml
```

This command:
- Computes MD5 checksums for all local files
- Compares against remote storage checksums
- Reports any mismatches (indicates corrupted files)

### Policy Management Commands

#### Generate Access Policy

Generate an IAM policy for benchmark access:

```bash
ob remote policy create -b benchmark.yaml
```

Example output:
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": ["s3:*"],
      "Resource": [
        "arn:aws:s3:::benchmark-name/*",
        "arn:aws:s3:::benchmark-name"
      ]
    },
    {
      "Effect": "Deny",
      "Action": [
        "s3:BypassGovernanceRetention",
        "s3:DeleteObjectTagging",
        "s3:DeleteObjectVersion",
        "s3:DeleteObjectVersionTagging",
        "s3:DeleteBucket",
        "s3:DeleteBucketPolicy"
      ],
      "Resource": [
        "arn:aws:s3:::benchmark-name/*",
        "arn:aws:s3:::benchmark-name"
      ]
    }
  ]
}
```

Use this output when creating access keys (see [Creating Access Keys](#creating-access-keys)).

### Archive Commands

#### Create Archive

Create a complete benchmark archive:

```bash
# Archive everything (config, code, software, results)
ob archive -b benchmark.yaml --code --software --results

# With compression
ob archive -b benchmark.yaml --code --software --results --compression lzma

# Dry run (preview what will be archived)
ob archive -b benchmark.yaml --code --software --results --dry-run
```

Archive components can be selected individually:
- `--code`: Include Git repositories
- `--software`: Include software environment definitions
- `--results`: Include benchmark output files

Compression levels:
- `--compression none`: No compression (fastest)
- `--compression deflated`: Standard compression (default)
- `--compression bzip2`: Better compression
- `--compression lzma`: Best compression (slowest)

The archive will be created in the current directory with a name like:
```
benchmark-name-0.1.zip
```

## Configuration

### Environment Variables

Storage credentials can be configured via environment variables:

```bash
# Option 1: Config file path
export OB_STORAGE_S3_CONFIG=/path/to/config.json

# Option 2: Direct credentials
export OB_STORAGE_S3_ACCESS_KEY=your_access_key
export OB_STORAGE_S3_SECRET_KEY=your_secret_key

# Benchmark path (optional, to avoid -b flag)
export OB_BENCHMARK=/path/to/benchmark.yaml
```

### Configuration File Format

Save credentials in JSON format:

```json
{
  "endpoint": "https://s3.example.com",
  "access_key": "YOUR_ACCESS_KEY",
  "secret_key": "YOUR_SECRET_KEY",
  "secure": true
}
```

**Note**: Keep this file secure and do not commit it to version control.

### Benchmark YAML Configuration

Specify storage endpoint in your benchmark YAML:

```yaml
id: my-benchmark
storage_api: S3
storage: https://s3.example.com
version: "0.1"
```

## Creating Access Keys

Access keys provide authentication for reading/writing benchmark results to remote storage.

### MinIO

Do note the latest minio community edition removed the access keys from the UI:
https://github.com/minio/minio/issues/21317#issuecomment-2910267585

You need to refer to the minio documentation for creating access keys:
https://min.io/docs/minio/linux/operations/create-access-key.html

### AWS

1. Create a new IAM user in AWS Console
2. Navigate to **IAM** → **Policies** → **Create policy**
3. Use the JSON editor and paste output from:
   ```bash
   ob remote policy create -b benchmark.yaml
   ```
4. Create the policy with a descriptive name
5. Attach the policy to the IAM user
6. Create an access key for the user
7. Save the access key and secret key

## Workflow Example

Complete workflow for using remote storage:

```bash
# 1. Generate access policy
ob remote policy create -b benchmark.yaml > policy.json

# 2. Create access key in MinIO/AWS using policy.json

# 3. Save credentials
echo '{"access_key": "KEY", "secret_key": "SECRET"}' > ~/.ob/s3-config.json

# 4. Run benchmark with remote storage
export OB_STORAGE_S3_CONFIG=~/.ob/s3-config.json
ob run benchmark -b benchmark.yaml

# 5. Create version snapshot
ob remote version create -b benchmark.yaml

# 6. List versions
ob remote version list -b benchmark.yaml

# 7. Download results
ob remote files download -b benchmark.yaml

# 8. Verify integrity
ob remote files checksum -b benchmark.yaml

# 9. Create archive
ob archive -b benchmark.yaml --code --software --results --compression lzma
```

## Related Documentation

- [Tutorial: Remote Storage](../../docs/src/tutorial.md#remote-storage---s3-aws-or-minio)
- [How-to Guides](../../docs/src/howto.md)
- [Benchmark YAML Reference](../../docs/src/bench-reference.md)

## See Also

- [MinIO Documentation](https://min.io/docs/minio/linux/index.html)
- [AWS S3 Documentation](https://docs.aws.amazon.com/s3/)
- [S3 Object Lock](https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lock.html)
