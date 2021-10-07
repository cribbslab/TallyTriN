__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from src.util.random.Ordering import ordering as ranord
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.sequencing.reads.simulate.pcr.Error import error as pcrerr


class amplify(object):

    def __init__(self, pcr_params):
        self.pcr_params = pcr_params

    def np(self, ):
        for ipcr in range(self.pcr_params['pcr_num']):
            print('--->at PCR {}'.format(ipcr + 1))
            self.pcr_params['ipcr'] = ipcr
            self.pcr_params = self.flow(params=self.pcr_params)
            # print(std_flow_params.keys())
        return self.pcr_params

    @pcrerr(method='default')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow(self, params):
        return params