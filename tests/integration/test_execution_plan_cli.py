#!/usr/bin/env python3
"""Integration tests for MetaNexus execution plan CLI invocations and side-effect guarantees."""

import unittest

try:
    import pytest
    pytestmark = pytest.mark.integration
except ImportError:
    pass

from tests.test_execution_plan import TestCLISideEffectsAndInvocations

if __name__ == "__main__":
    unittest.main()
