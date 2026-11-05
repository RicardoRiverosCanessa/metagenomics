try:
    import pytest
    pytestmark = pytest.mark.container
except ImportError:
    pass
