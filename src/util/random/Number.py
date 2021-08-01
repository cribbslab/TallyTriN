__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np
from functools import wraps


class number(object):

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, deal):
        if self.kwargs['type'] == 'binomial':
            distrib = self.binomial
        if self.kwargs['type'] == 'uniform':
            distrib = self.binomial
        else:
            pass
        @wraps(deal)
        def switch(ph, *args, **kwargs):
            print('numbering...')
            # print(kwargs)
            res = deal(ph, **kwargs)
            res['num'] = distrib(
                n=len(res['data']),
                p=res['p'],
            )
            # print(res)
            return res
        return switch

    def binomial(self, n, p, use_seed=True, seed=1):
        if use_seed:
            state = np.random.RandomState(seed)
            return state.binomial(
                n,
                p,
            )
        else:
            return np.random.binomial(
                n,
                p,
            )

    def uniform(self, low, high, num, use_seed=True, seed=1):
        if use_seed:
            state = np.random.RandomState(seed)
            return state.randint(
                low=low,
                high=high,
                size=num
            )
        else:
            return np.random.randint(
                low=low,
                high=high,
                size=num
            )