'''test_import - test importing all modules
===========================================

This script attempts to import all the python code
in the repository and tests for successful loading.

This script is best run within nosetests::

   pytest tests/test_import.py

'''
import importlib
import os
import pytest

# define the directories to test
directories = ['tallytrin']

# define the test function
@pytest.mark.parametrize('directory', directories)
def test_imports(directory):
    for file in directory:
        module = importlib.import_module(directory)
        assert module is not None
