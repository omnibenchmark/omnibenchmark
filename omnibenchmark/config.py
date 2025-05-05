"""Configuration to set up a local cache and a datadir for test data download"""

import os
import yaml

_home = os.path.expanduser("~")

app_name = "omni-py"

xdg_bench_home = os.environ.get("XDG_DATA_HOME") or os.path.join(
    _home, ".local", "share"
)

xdg_config_home = os.environ.get("XDG_CONFIG_HOME") or os.path.join(_home, ".config")

bench_dir = os.path.join(xdg_bench_home, app_name)
config_dir = os.path.join(xdg_config_home, app_name)

rc_file = os.path.join(config_dir, "omni-py.yaml")

default_cfg = {"dirs": {"datasets": "~/OmniBenchmark/datasets"}}


def init_dirs():
    os.makedirs(bench_dir, exist_ok=True)
    os.makedirs(config_dir, exist_ok=True)


def init_rc():
    if not os.path.isfile(rc_file):
        _write_config(default_cfg)


def get_dataset_dir():
    c = _get_config()
    path = c.get("dirs").get("datasets")
    return os.path.expanduser(path)


def _get_config():
    with open(rc_file) as f:
        return yaml.safe_load(f.read())


def _write_config(c):
    with open(rc_file, "w") as f:
        yaml.dump(c, f)
