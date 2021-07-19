__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd


class group(object):

    def __init__(self, ):
        pass

    def uniform(self, df, num, use_seed=True, seed=1, replace=False):
        df_ = pd.DataFrame(df)
        num_samples = df_.shape[0]
        if use_seed:
            state = np.random.RandomState(seed)
            df_ = df_.iloc[state.choice(num_samples, num, replace=replace)]
        else:
            df_ = df_.iloc[np.random.choice(num_samples, num, replace=replace)]
        return df_