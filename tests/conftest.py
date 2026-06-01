import pytest
import os

@pytest.fixture(scope="session")
def data_dir():
    return os.path.join(os.path.dirname(__file__), "data")
