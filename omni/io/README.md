
# How to use S3Storage/MinIOStorage

## get new access key

In the MinIO Console navigate to 'Access Keys' and click 'Create access key'. Set 'Restrict beyond user policy' to 'ON'. Create new policy with

```
import omni
import json
policy = omni.io.S3config.benchmarker_access_token_policy(<BENCHMARK>)
print(json.dumps(policy, indent=2))
```

Replace the displayed policy with the output of the above command.
Enter a name (Preferable syntax: `Benchmark: <BENCHMARK>`) and a description. Click on `Create` and copy the access key and secret key. Save in a `<CONFIG>.json` file somewhere with the following format:

```
{"endpoint": "<URL>", "access_key": "<ACCESS_KEY>", "secret_key": "<SECURE_KEY>", "secure": false}
```

## usage

The to use the MinIOStorage or S3Storage class, use the following code:

```
import json
from omni.io.MinIOStorage import MinIOStorage

with open("<CONFIG>.json", "r") as file:
    auth_options = json.load(file)

ss = MinIOStorage(auth_options=auth_options, benchmark="<BENCHMARK>")
ss.versions
```
