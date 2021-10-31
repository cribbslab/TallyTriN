__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from simreadflow.sequencing.Error import error as seqerr
from simreadflow.util.random.Sampling import sampling as ranspl


class calling(object):

    def __init__(self, seq_params):
        self.seq_params = seq_params

    def np(self, ):
        self.seq_params['spl_num'] = int(self.seq_params['data'].shape[0] * self.seq_params['seq_sub_spl_rate'])
        res_flow = self.flow(params=self.seq_params)
        # print(res_flow)
        return res_flow

    @seqerr(method='default')
    @ranspl(method='uniform')
    def flow(self, params):
        return params