__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from functools import wraps


class library(object):

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, deal):
        manage = self.kwargs['method']
        @wraps(deal)
        def build(ph, *args, **kwargs):
            res = deal(ph, **kwargs)
            if kwargs['is_sv'] is True:
                if manage == 'default':
                    with open(kwargs['lib_fpn'], 'a') as file:
                        file.write(res + "\n")
                elif manage == 'separate':
                    with open(kwargs['lib_fpn'], 'a') as file:
                        file.write(kwargs['res'] + "\n")
            return res
        return build