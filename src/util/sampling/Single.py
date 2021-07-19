__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np


class single(object):

    def __init__(self, ):
        pass

    def uniform(self, low, high, num):
        return np.random.randint(low=low, high=high, size=num)


if __name__ == "__main__":
    p = single()
    print(p.uniform(
        low=0,
        high=4,
        num=60,
    ))