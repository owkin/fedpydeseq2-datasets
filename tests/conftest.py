import json
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def raw_data_path():
    """Fixture to get the path to the raw data."""
    default_paths = Path(__file__).parent / "paths_default.json"
    specified_paths = Path(__file__).parent / "paths.json"
    if specified_paths.exists():
        with open(specified_paths) as f:
            raw_data_path = json.load(f)["raw_data"]
    else:
        with open(default_paths) as f:
            raw_data_path = json.load(f)["raw_data"]
    if raw_data_path.startswith("/"):
        raw_data_path = Path(raw_data_path)
    else:
        raw_data_path = Path(__file__).parent / raw_data_path
    print("Test raw data path")
    return raw_data_path


@pytest.fixture(scope="session")
def tmp_processed_data_path(tmpdir_factory):
    return Path(tmpdir_factory.mktemp("processed"))
