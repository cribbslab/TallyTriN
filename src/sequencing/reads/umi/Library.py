__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np
from functools import wraps


class library(object):

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, deal):
        @wraps(deal)
        def build(ph, *args, **kwargs):
            res = deal(ph, **kwargs)
            with open(kwargs['umi_lib_fpn'], 'a') as file:
                file.write(res + "\n")
            return res
        return build