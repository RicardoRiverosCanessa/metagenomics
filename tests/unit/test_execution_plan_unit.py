#!/usr/bin/env python3
"""Unit tests for MetaNexus execution plan model, schema, graph, and invariants."""

import unittest

try:
    import pytest
    pytestmark = pytest.mark.unit
except ImportError:
    pass

from tests.test_execution_plan import (
    TestApplicabilityAndPrecedence,
    TestCheckpointsAndValidation,
    TestExecutabilityAndSummary,
    TestForceAndRestart,
    TestGraphAndInvalidation,
    TestModelSchemaValidation,
)

if __name__ == "__main__":
    unittest.main()
