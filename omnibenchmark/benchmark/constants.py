import os
from typing import Final

# LOCAL_TIMEOUT_VAR is the name of the environment variable that will be used to set the timeout
# for the spawned task on local environments. Snakemake will act weirdly if we reuse either 'runtime' or 'timeout',
# so don't be tempted to change that.
LOCAL_TIMEOUT_VAR: Final[str] = "local_task_timeout"

OUTPUT_PATH_PREFIX: Final[str] = os.path.join(
    "{input}", "{stage}", "{module}", "{params}"
)
