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
directories = ['tallytrin/entry.py']

# define the test function
@pytest.mark.parametrize('directory', directories)
def test_imports(directory):
    for file in os.listdir(directory):
        if file.endswith('.py') and file != '__init__.py':
            module_name = file[:-3]
            module = importlib.import_module(f'{directory}.{module_name}')
            assert module is not None
