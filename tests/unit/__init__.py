try:
    import pytest
    pytestmark = pytest.mark.unit
except ImportError:
    pass
