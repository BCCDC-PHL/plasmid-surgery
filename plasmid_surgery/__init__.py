"""Assess plasmid reconstructions against reference plasmids."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("plasmid-surgery")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0+unknown"
