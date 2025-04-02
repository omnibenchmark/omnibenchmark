from typing import Final

# DEFAULT_TIMEOUT_SECONDS is the max value we'll try to enforce as a timeout for the spawned task.
# This is currently only enforced on local environments (but we can pass the same value to job queues)
DEFAULT_TIMEOUT_SECONDS: Final[int] = 8 * 60 * 60  # Default is 4 hours
