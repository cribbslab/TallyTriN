import pycodestyle
import glob
import os

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


def check_style(filename):
    '''check style of filename.'''

    style_guide = pycodestyle.StyleGuide(quiet=True, ignore=IGNORE)
    report = style_guide.check_files([filename])

    # count errors/warnings excluding those to ignore
    assert report.total_errors == 0, f"Style violations in {filename}"


def test_style():
    '''test style of scripts'''

    for label, expression in EXPRESSIONS:
        files = glob.glob(expression)
        files.sort()

        for f in files:
            if os.path.isdir(f):
                continue
            yield check_style, os.path.abspath(f)
