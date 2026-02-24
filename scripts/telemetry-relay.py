#!/usr/bin/env python3
"""
Thin shim — delegates to omnibenchmark.telemetry.relay.

Prefer running via the package:
    python -m omnibenchmark.telemetry.relay [options]
"""

from omnibenchmark.telemetry.relay import main

if __name__ == "__main__":
    main()
