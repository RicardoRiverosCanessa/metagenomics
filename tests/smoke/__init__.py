try:
    import pytest
    pytestmark = pytest.mark.smoke
except ImportError:
    pass
