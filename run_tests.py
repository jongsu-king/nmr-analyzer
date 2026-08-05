#!/usr/bin/env python3
"""Run the whole test suite: python3 run_tests.py"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))

if __name__ == "__main__":
    suite = unittest.defaultTestLoader.discover(os.path.join(HERE, "tests"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
