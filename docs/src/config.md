# Configuration System

Omnibenchmark uses a configuration file to store settings such as paths and environment variables. This document describes both the user-facing aspects and the programmatic interface for developers.

## File Location

The configuration file is stored at:

- Linux: `~/.config/omnibenchmark/omnibenchmark.cfg`
- macOS: `~/Library/Application Support/omnibenchmark/omnibenchmark.cfg`

## Format

The configuration file uses the INI format with sections and key-value pairs:

```ini
[section1]
key1 = value1
key2 = value2

[section2]
key3 = value3
```

## Common Sections and Keys

Omnibenchmark uses the following standard configuration sections:

- `easybuild`: Contains settings for Easybuild integration
  - `MODULEPATH`: Path to module installations
  - `ROBOTPATH`: Path to easyconfigs repository

- `dirs`: Contains paths for dataset storage and git module caching
  - `datasets`: Path to store benchmark datasets
  - `git_modules`: Path to cache cloned git repositories (default: `.omnibenchmark/git`)

## Programmatic Access

The configuration system provides a Python API for reading and writing configuration values.

### Reading Configuration Values

```python
from omnibenchmark.config import config

# Get a value with a default if key doesn't exist
value = config.get('section', 'key', default='default_value')
```

### Writing Configuration Values

```python
from omnibenchmark.config import config

# Set a value
config.set('section', 'key', 'new_value')

# Save changes to disk
config.save()
```

### Accessing Configuration Options

```python
from omnibenchmark.config import config

# Get all sections
sections = config.sections()

# Get all keys in a section
keys = config.options('section')
```

### Creating a Custom Configuration Accessor

In most cases, you should use the global `config` instance. However, you can create a custom accessor if needed:

```python
from omnibenchmark.config import ConfigAccessor
from pathlib import Path

# Create a custom config accessor with a specific path
custom_config = ConfigAccessor(Path('/path/to/custom/omnibenchmark.cfg'))

# Use it the same way as the global config
value = custom_config.get('section', 'key', default='default_value')
```

## Implementation Details

The configuration system is implemented in `omnibenchmark/config.py`. It provides:

- Platform-specific paths for configuration files
- A `ConfigAccessor` class that provides dict-like access to the configuration
- Initialization of necessary directories
- Graceful handling of missing sections or keys

When adding new configuration options to Omnibenchmark, be sure to update this documentation with the new sections and keys.