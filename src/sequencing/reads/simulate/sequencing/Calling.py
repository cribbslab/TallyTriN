__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.sequencing.reads.simulate.sequencing.Error import error as seqerr


class calling(object):

    def __init__(self, calling_params):
        self.calling_params = calling_params

    def np(self, ):
        std_flow_params = self.flow(params=self.calling_params)
        return std_flow_params

    @seqerr(method='tomap')
    def flow(self, params):
        return params