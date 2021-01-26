from bodas import __version__
import toml

def _read_version_from_toml_file():
    with open("pyproject.toml", "r") as f:
        content = f.read()
    return toml.loads(content)['tool']['poetry']['version']


def test_version():
    assert __version__ == _read_version_from_toml_file()
