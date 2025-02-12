import io

import pytest
import logging


@pytest.fixture
def capture_logs():
    """Fixture to capture log output during tests."""
    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)
    logger = logging.getLogger("omnibenchmark")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    yield log_stream  # Yield the stream to the test function

    # Cleanup after the test
    logger.removeHandler(handler)
    log_stream.close()
