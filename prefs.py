"""Persistent application preferences: recent files, window geometry, defaults.

Stored as JSON in the usual per-platform location.  Every read is defensive —
a corrupt or hand-edited preferences file must never stop the program from
starting.
"""

from __future__ import annotations

import json
import os
import sys

APP_DIR_NAME = "NMR Analyzer"
MAX_RECENT = 10

DEFAULTS = {
    "recent": [],
    "last_dir": "",
    "geometry": "",
    "sensitivity": 8.0,
    "line_broadening": 0.3,
    "zero_fill": "65536",
    "baseline_method": "spline",
    "stack": False,
    "normalise_each": True,
    "show_peaks": True,
    "show_integrals": True,
    "show_grid": False,
}


def config_dir():
    if sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Application Support")
    elif os.name == "nt":
        base = os.environ.get("APPDATA", os.path.expanduser("~"))
    else:
        base = os.environ.get("XDG_CONFIG_HOME",
                              os.path.expanduser("~/.config"))
    return os.path.join(base, APP_DIR_NAME)


def config_path():
    return os.path.join(config_dir(), "preferences.json")


class Preferences:
    def __init__(self):
        self._data = dict(DEFAULTS)
        self.load()

    # -- persistence --------------------------------------------------------

    def load(self):
        try:
            with open(config_path(), "r", encoding="utf-8") as fh:
                stored = json.load(fh)
            if isinstance(stored, dict):
                for key, value in stored.items():
                    if key in DEFAULTS:
                        self._data[key] = value
        except (OSError, ValueError):
            pass                      # first run, or the file is unreadable
        if not isinstance(self._data.get("recent"), list):
            self._data["recent"] = []

    def save(self):
        try:
            os.makedirs(config_dir(), exist_ok=True)
            with open(config_path(), "w", encoding="utf-8") as fh:
                json.dump(self._data, fh, indent=2, ensure_ascii=False)
        except OSError:
            pass                      # preferences are a convenience, not data

    # -- access -------------------------------------------------------------

    def get(self, key, default=None):
        return self._data.get(key, DEFAULTS.get(key, default))

    def set(self, key, value):
        self._data[key] = value

    # -- recent files -------------------------------------------------------

    @property
    def recent(self):
        """Recently used paths, most recent first, existing ones only."""
        return [p for p in self._data.get("recent", []) if os.path.exists(p)]

    def add_recent(self, path):
        path = os.path.abspath(path)
        items = [p for p in self._data.get("recent", [])
                 if os.path.abspath(p) != path]
        items.insert(0, path)
        self._data["recent"] = items[:MAX_RECENT]
        self._data["last_dir"] = os.path.dirname(path)

    def clear_recent(self):
        self._data["recent"] = []
