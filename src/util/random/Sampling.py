__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
from functools import wraps


class sampling(object):

    def __init__(self, method='uniform'):
        self.method = method

    def __call__(self, deal):
        if self.method == 'uniform':
            compiler = self.uniform
        elif self.method == 'gaussian':
            compiler = self.gaussian
        elif self.method == 'poisson':
            compiler = self.poisson
        else:
            compiler = self.uniform
        @wraps(deal)
        def switch(dself, *args, **kwargs):
            # print(args)
            # print(kwargs)
            print('------>sampling...')
            res2p = deal(dself, **kwargs)
            # print(res2p)
            res2p['data_spl'] = compiler(
                data=res2p['data'],
                num=res2p['ipcr_num'],
                # use_seed=res2p['use_seed'],
                # seed=res2p['seed'],
                # replace=res2p['replace'],
            )
            # print(res2p['data'])
            return res2p
        return switch

    def uniform(self, data, num, use_seed=True, seed=1, replace=False):
        """
        uniform
        :param data:
        :param num:
        :param use_seed:
        :param seed:
        :param replace:
        :return:
        """
        num_samples = len(data)
        if isinstance(data, pd.DataFrame):
            if use_seed:
                state = np.random.RandomState(seed)
                data = data.iloc[state.choice(num_samples, num, replace=replace)].reset_index(drop=True)
            else:
                data = data.iloc[np.random.choice(num_samples, num, replace=replace)].reset_index(drop=True)
        elif type(data) is np.ndarray:
            if use_seed:
                state = np.random.RandomState(seed)
                data = data[state.choice(num_samples, num, replace=replace)]
            else:
                data = data[np.random.choice(num_samples, num, replace=replace)]
        else:
            if use_seed:
                state = np.random.RandomState(seed)
                ids = state.choice(num_samples, num, replace=replace)
                data = [data[i - 1] for i in ids]
            else:
                ids = np.random.choice(num_samples, num, replace=replace)
                data = [data[i - 1] for i in ids]
        return data

    def gaussian(self, ):
        pass
    
    def poisson(self, ):
        pass