__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
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
            print('======>numbering...')
            # print(kwargs)
            res = deal(ph, **kwargs)
            res['spl_num'] = distrib(
                n=len(res['data']),
                p=res['ampl_rate'],
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

    def nbinomial(self, n, p, use_seed=True, seed=1):
        """

        :param n: the number of success to be expected, better, n = the total number of trails * p
        :param p: the prob of success
        :param use_seed:
        :param seed:
        :return:
        """
        if use_seed:
            state = np.random.RandomState(seed)
            return state.negative_binomial(
                n,
                p,
            )
        else:
            return np.random.negative_binomial(
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