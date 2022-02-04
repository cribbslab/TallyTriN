__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
from functools import wraps


class roller(object):

    def __init__(self, method):
        self.method = method

    def __call__(self, deal):
        if self.method == 'uniform':
            sample = self.uniform
        elif self.method == 'gaussian':
            sample = self.gaussian
        else:
            sample = self.uniform

        @wraps(deal)
        def switch(dself, *args, **kwargs):
            # print(args)
            # print(kwargs)
            data = deal(dself, **kwargs)
            return sample(
                data=data,
                num=kwargs['data'],
                use_seed=kwargs['data'],
                seed=kwargs['data'],
                replace=kwargs['data'],
            )

        return switch

    def uniform(self, data, num, use_seed=True, seed=1, replace=False):
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