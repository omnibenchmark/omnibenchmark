from typing import Final

from humanfriendly import parse_timespan

# DEFAULT_TIMEOUT_HUMAN is the human representation of the max value we'll try to enforce
# as a timeout for the spawned task.
# This is currently only enforced on local environments; when supporting slurm we can pass
# it as runtime resource.
DEFAULT_TIMEOUT_HUMAN: Final[str] = "4h"
DEFAULT_TIMEOUT_SECONDS: Final[int] = int(parse_timespan(DEFAULT_TIMEOUT_HUMAN))

# LOCAL_TIMEOUT_VAR is the name of the environment variable that will be used to set the timeout
# for the spawned task on local environments. Snakemake will act weirdly if we reuse either 'runtime' or 'timeout',
# so don't be tempted to change that.
LOCAL_TIMEOUT_VAR: Final[str] = "local_task_timeout"
