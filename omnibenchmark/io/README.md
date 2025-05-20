# Overview

Helper functions to interact with remote storage (default: MinIO) to create versions of benchmarks. The following is a short description.
For each benchmark a version-controlled S3 bucket is created. Files from a benchmark run can be uploaded to this bucket. For example:
```
/
|-- out
    |-- f1.txt
    |-- f2.txt
```
The versioning means that new uploads of files with existing names get a new versionID. For example:
```
/
|-- out
    |-- f1.txt
        (versionID: V2)
        (versionID: V1)
    |-- f2.txt
        (versionID: V2)
        (versionID: V1)
```
To create a benchmark version (or snapshot) the newest file versions in specified directories (default: 'out','versions') get tagged with the benchmark version and write protected (setting retention policy to Governance on single file versions). Additionally an overview file gets created resulting in the following structure for tagging with `version=0.1`:
```
/
|-- versions
    |-- 0.1.csv
|-- out
    |-- f1.txt
        (versionID: V2, tags: 0.1, Retention Policy: Governance)
        (versionID: V1)
    |-- f2.txt
        (versionID: V2, tags: 0.1, Retention Policy: Governance)
        (versionID: V1)
```
The `0.1.csv` file contains columns `name`,`version_id`,`last_modified`,`size`,`etag`.


# How to use MinIOStorage (S3)

## Create policy

Create new policy with
```
import omnibenchmark
import json
policy = omnibenchmark.io.S3config.benchmarker_access_token_policy(<BENCHMARK>)
print(json.dumps(policy, indent=2))
```

## Create new access key

### MinIO

In the MinIO Console navigate to 'Access Keys' and click 'Create access key'. Set 'Restrict beyond user policy' to 'ON'. Replace the displayed policy with the output of the above command.
Enter a name (Preferable syntax: `Benchmark: <BENCHMARK>`) and a description. Click on `Create` and copy the access key and secret key.

### AWS

Create a new user. Create a new policy with the output of the above command. Attach policy to user. Create access key for user.

## Save access key information locally

Save the access key and secret key in a `<CONFIG>.json` file somewhere with the following format:

```
{"endpoint": "<URL>", "access_key": "<ACCESS_KEY>", "secret_key": "<SECURE_KEY>", "secure": false}
```

## Usage

To use the MinIOStorage class, use the following code:

```
import json
from omnibenchmark.io.MinIOStorage import MinIOStorage

with open("<CONFIG>.json", "r") as file:
    auth_options = json.load(file)

ms = MinIOStorage(auth_options=auth_options, benchmark="<BENCHMARK>")
ms.versions
```

## Example usage

See [examples/storage_usage.py](examples/storage_usage.py).
