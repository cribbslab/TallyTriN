import pycodestyle
import glob
import os
import pytest

# DIRECTORIES to examine
EXPRESSIONS = (
    ('FirstLevel', 'scpipelines/*.py'),
    ('SecondLevel', 'scpipelines/version.py'))

# Codes to ignore in the pycodestyle BaseReport
IGNORE = {
    'E101',  # indentation contains mixed spaces and tabs
    'E201',  # whitespace after '('
    'E202',  # whitespace before ')'
    'E122',  # continuation line missing indentation or outdented
    'E265',  # block comment should start with '# '
    'E501',  # line too long (82 > 79 characters)
    'E502',  # the backslash is redundant between brackets
    'E731',  # do not assign a lambda expression, use a def
    'W191',
    'W291',
    'W293',
    'W391',
    'W503',  # line break before binary operator
    'W601',
    'W602'
}


def _collect_style_files():
    '''collect all Python files to check.'''
    files = []
    for label, expression in EXPRESSIONS:
        for f in sorted(glob.glob(expression)):
            if not os.path.isdir(f):
                files.append(os.path.abspath(f))
    return files


@pytest.mark.parametrize("filename", _collect_style_files())
def test_style(filename):
    '''test style of a script.'''

    style_guide = pycodestyle.StyleGuide(quiet=True, ignore=IGNORE)
    report = style_guide.check_files([filename])

    assert report.total_errors == 0, f"Style violations in {filename}"
