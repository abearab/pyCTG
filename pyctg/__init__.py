from .main import (
    CTG_synergy, 
    read_CTG_titration_data, 
    read_CTG_synergy_data
)


def _get_version():

    import os

    pyproject_path = os.path.join(os.path.dirname(__file__), "..", "pyproject.toml")

    with open(pyproject_path, "r") as pyproject_file:
        for line in pyproject_file.readlines():
            if "version" in line:
                return line.split("=")[1].strip().strip('"')


try:
    __version__ = _get_version()
except Exception:
    __version__ = "Unknown"
