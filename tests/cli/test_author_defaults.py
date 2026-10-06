import subprocess

import pytest

from omnibenchmark.cli import create
from omnibenchmark.config import ConfigAccessor


@pytest.mark.short
def test_author_defaults_precedence(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    monkeypatch.chdir(tmp_path)
    cfg = ConfigAccessor(tmp_path / "ob.cfg")
    monkeypatch.setattr(create, "config", cfg)

    assert create._author_defaults() == {}

    subprocess.run(["git", "init", "-q"], check=True)
    subprocess.run(["git", "config", "user.name", "Git Name"], check=True)
    subprocess.run(["git", "config", "user.email", "git@x.org"], check=True)
    assert create._author_defaults() == {
        "author_name": "Git Name",
        "author_email": "git@x.org",
    }

    cfg.set("user", "name", "Ob Name")
    assert create._author_defaults() == {
        "author_name": "Ob Name",
        "author_email": "git@x.org",
    }

    assert cfg.unset("user", "name")
    assert not cfg.unset("user", "name")
    assert "user" not in cfg.sections()
    assert create._author_defaults()["author_name"] == "Git Name"
