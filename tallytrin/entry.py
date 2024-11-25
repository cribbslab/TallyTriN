'''
TallyTriN.py - collection of bulk and single-cell workflows 
===========================================================

:Tags: Single-cell

To use a specific workflow, type::

    tallytrin <workflow> [workflow options] [workflow arguments]

For this message and a list of available keywords type::

    tallytrin --help

To get help for a specific workflow, type::

    tallytrin <workflow> --help
'''

import os
import sys
import re
import glob
import importlib.util
import tallytrin


def print_list_in_columns(l, ncolumns):
    '''Output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % ncolumns != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):

    argv = sys.argv

    # paths to look for pipelines:
    path = os.path.abspath(os.path.dirname(tallytrin.__file__))
    relpath = os.path.abspath("../src")

    paths = [path, relpath]

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        pipelines = []
        for path in paths:
            pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of available pipelines are:\n")
        print("{}\n".format(
            print_list_in_columns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                2)))
        return

    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)

    # remove 'tallytrin' from sys.argv
    del sys.argv[0]

    spec = None
    for path in paths:
        try:
            spec = importlib.util.spec_from_file_location(pipeline, os.path.join(path, f"{pipeline}.py"))
            if spec is not None:
                break
        except FileNotFoundError:
            continue

    if spec is None or spec.loader is None:
        print(f"Error: pipeline '{command}' not found.")
        sys.exit(1)

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
