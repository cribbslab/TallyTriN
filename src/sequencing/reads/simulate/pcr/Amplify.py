__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.util.random.Ordering import ordering as ranord
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.sequencing.reads.simulate.pcr.Error import error as pcrerr


class amplify(object):
    
    def __init__(self, pcr_params):
        self.pcr_params = pcr_params

    def np(self, ):
        import numpy as np
        std_flow_params = {
            'data': np.array(self.pcr_params['init_seqs']),
            'p': self.pcr_params['ampl_rate'],
            'pcr_error': self.pcr_params['pcr_error'],
            'use_seed': self.pcr_params['use_seed'],
            'seed': self.pcr_params['seed'],
            'recorder_num_bases': [],
            'recorder_num_pcr_errs': [],
            'recorder_num_pcr_reads': [],
        }
        for ipcr in range(self.pcr_params['num_pcr']):
            std_flow_params['ipcr'] = ipcr
            std_flow_params = self.flow(params=std_flow_params)
            # print(std_flow_params.keys())
        return std_flow_params

    @pcrerr(method='tomap')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow(self, params):
        return params

    # @errroller(method='uniform')
    # def splUniform(self, data):
    #     return data