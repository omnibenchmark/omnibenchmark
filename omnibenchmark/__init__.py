try:
    from importlib.metadata import version

    __version__ = version("omnibenchmark")

except Exception:
    __version__ = "0.0.0"  # Fallback version for development
