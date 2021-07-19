__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

import sys
sys.path.append('../../../../../')


class single(object):

    def __init__(self, ):
        pass

    def get(self, gap=False, universal=False):
        if universal:
            if gap:
                return ['A', 'C', 'G', 'T', '-']
            else:
                return ['A', 'C', 'G', 'T']
        else:
            if gap:
                return ['A', 'T', 'C', 'G', '-']
            else:
                return ['A', 'T', 'C', 'G']

    def todict(self, bases, reverse=False):
        aa_dict = {}
        for k, v in enumerate(bases):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict


if __name__ == "__main__":
    p = single()
    bs = p.get()
    print(bs)
    print(p.todict(bs, reverse=True))