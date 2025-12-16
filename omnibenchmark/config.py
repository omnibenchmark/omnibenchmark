"""Configuration to set up a local cache and a datadir for test data download"""

# TODO: make sure this is in use

import configparser
import os
import platform
from typing import Optional, Any

from pathlib import Path

APP_NAME = "omnibenchmark"

_home = os.path.expanduser("~")

xdg_config_home = os.environ.get("XDG_CONFIG_HOME") or os.path.join(_home, ".config")
xdg_bench_home = os.environ.get("XDG_DATA_HOME") or os.path.join(
    _home, ".local", "share"
)


default_cfg = {
    "dirs": {"datasets": f"~/{APP_NAME}/datasets", "git_modules": f".{APP_NAME}/git"}
}

bench_dir = os.path.join(xdg_bench_home, APP_NAME)

if platform.system() == "Darwin":
    # macOS
    config_dir = Path("~/Library/Application Support/omnibenchmark").expanduser()
else:
    # Linux or others
    config_dir = Path(os.path.join(xdg_config_home, APP_NAME))


def get_config_file():
    return config_dir / f"{APP_NAME}.cfg"


def init_dirs():
    """Initialize configuration and benchmark directories.

    Fails gracefully if directories cannot be created (e.g., read-only filesystem).
    """
    import logging

    logger = logging.getLogger(__name__)

    try:
        os.makedirs(bench_dir, exist_ok=True)
    except OSError as e:
        logger.warning(
            f"Could not create benchmark directory {bench_dir}: {e}. "
            "This may affect git module caching but won't prevent execution."
        )

    try:
        os.makedirs(config_dir, exist_ok=True)
    except OSError as e:
        logger.warning(
            f"Could not create config directory {config_dir}: {e}. "
            "Using in-memory configuration only."
        )


# def _get_config(config_path):
#    config = configparser.ConfigParser()
#    config.read(config_path)
#    return config


class ConfigAccessor:
    """
    A dict-like accessor for configuration files.

    This class provides a way to access configuration options with a dictionary-like
    interface while handling missing sections or keys gracefully.

    Usage:
        config = ConfigAccessor()
        value = config.get('section', 'key', fallback='default')
    """

    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize a ConfigAccessor with an optional config file path.

        Args:
            config_path: Path to the configuration file. If None, uses the default path.
        """
        if config_path is None:
            self.config_path = get_config_file()
        else:
            self.config_path = config_path

        # Ensure the config directory exists
        init_dirs()

        # Initialize the config parser
        self.config = configparser.ConfigParser()
        if self.config_path.exists():
            self.config.read(self.config_path)

    def get(self, section: str, key: str, default: Any = None) -> Any:
        """
        Get a configuration value from the specified section and key.

        Args:
            section: The configuration section
            key: The configuration key
            default: Value to return if the section or key doesn't exist

        Returns:
            The configuration value if it exists, otherwise the default value
        """
        try:
            return self.config[section][key]
        except (KeyError, configparser.NoSectionError, configparser.NoOptionError):
            return default

    def set(self, section: str, key: str, value: str) -> None:
        """
        Set a configuration value.

        Args:
            section: The configuration section
            key: The configuration key
            value: The value to set
        """
        if not self.config.has_section(section):
            self.config.add_section(section)

        self.config[section][key] = value

    def save(self) -> None:
        """
        Save the current configuration to the config file.

        Fails gracefully if the file cannot be written (e.g., read-only filesystem).
        """
        import logging

        logger = logging.getLogger(__name__)

        try:
            with open(self.config_path, "w") as configfile:
                self.config.write(configfile)
        except (OSError, IOError) as e:
            logger.warning(
                f"Could not save configuration to {self.config_path}: {e}. "
                "Configuration changes will not persist."
            )

    def sections(self) -> list:
        """
        Get all available sections in the config.

        Returns:
            List of section names
        """
        return self.config.sections()

    def options(self, section: str) -> list:
        """
        Get all options (keys) in a section.

        Args:
            section: The section name

        Returns:
            List of options in the section or empty list if section doesn't exist
        """
        try:
            return self.config.options(section)
        except configparser.NoSectionError:
            return []


# Create a global config accessor instance
config = ConfigAccessor()


def get_git_modules_dir() -> Path:
    """
    Get the configured directory for caching git modules/repositories.

    Returns:
        Path to the git modules cache directory (defaults to .omnibenchmark/git)
    """
    git_modules_dir_str = config.get(
        "dirs", "git_modules", default_cfg["dirs"]["git_modules"]
    )
    git_modules_dir = Path(git_modules_dir_str).expanduser()

    # Create directory if it doesn't exist
    git_modules_dir.mkdir(parents=True, exist_ok=True)

    return git_modules_dir
