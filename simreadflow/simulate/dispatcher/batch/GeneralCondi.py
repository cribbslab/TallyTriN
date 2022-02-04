__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from simreadflow.simulate.dispatcher.single.GeneralCondi import generalCondi as simugcondi
from Path import to


class generalCondi(object):

    def __init__(self, ):
        self.permutation_num = 10

        self.umi_unit_len_fixed = 10
        self.spacer_len_fixed = 10
        self.umi_num_fixed = 50
        self.pcr_num_fixed = 10
        self.pcr_err_fixed = 1e-3
        self.seq_err_fixed = 1e-3
        self.ampl_rate_fixed = 0.85

        self.ampl_rates = np.linspace(0.1, 1, 10)
        self.umi_unit_lens = np.arange(6, 36 + 1, 1)
        self.umi_nums = np.arange(20, 140 + 20, 20)
        self.pcr_nums = np.arange(1, 20 + 1, 1)
        self.pcr_errs, self.seq_errs = self.errors()
        print(self.pcr_errs, self.seq_errs)

        self.metrics = {
            'pcr_nums': self.pcr_nums,
            'pcr_errs': self.pcr_errs,
            'seq_errs': self.seq_errs,
            'ampl_rates': self.ampl_rates,
            'umi_lens': self.umi_unit_lens,
        }
        self.fastq_fn_pref = {
            'pcr_nums': 'pcr_',
            'pcr_errs': 'pcr_err_',
            'seq_errs': 'seq_err_',
            'ampl_rates': 'ampl_rate_',
            'umi_lens': 'umi_len_',
        }

    def errors(self, ):
        pcr_errs = []
        seq_errs = []
        e = 1e-5
        while e < 3e-1:
            pcr_errs.append(e)
            seq_errs.append(e)
            if 5 * e < 3e-1:
                pcr_errs.append(2.5 * e)
                pcr_errs.append(5 * e)
                pcr_errs.append(7.5 * e)
                seq_errs.append(2.5 * e)
                seq_errs.append(5 * e)

                seq_errs.append(7.5 * e)
            e = 10 * e
        pcr_errs.append(0.2)
        seq_errs.append(0.2)
        pcr_errs.append(0.3)
        seq_errs.append(0.3)
        # print(pcr_errs)
        # print(seq_errs)
        return pcr_errs, seq_errs

    def umiLens(self, ):
        for id, umi_len in enumerate(self.umi_unit_lens):
            simu_params = {
                'init_seq_setting': {
                    'seq_num': self.umi_num_fixed,
                    'umi_unit_pattern': 1,
                    'umi_unit_len': umi_len,
                    'is_seed': True,
                    'seq_len': 100 - umi_len - self.spacer_len_fixed,
                    'spacer_len': self.spacer_len_fixed,
                    'is_sv_umi_lib': True,
                    'is_sv_seq_lib': True,
                    'is_sv_spacer_lib': True,
                    'working_dir': to('data/simu/umi_lens/spacer/'),
                    'umi_lib_fpn': to('data/simu/umi_lens/spacer/umi_') + str(umi_len) + '.txt',
                    'seq_lib_fpn': to('data/simu/umi_lens/spacer/seq_') + str(umi_len) + '.txt',
                    'spacer_lib_fpn': to('data/simu/umi_lens/spacer/spacer_') + str(umi_len) + '.txt',
                    'condis': ['umi', 'spacer', 'seq'],
                    'permutation': 0
                },
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': self.pcr_num_fixed,
                'pcr_error': self.pcr_err_fixed,
                'seq_error': self.seq_err_fixed,
                'err_num_met': 'nbinomial',
                'use_seed': False,
                'seed': None,
                'write': {
                    'fastq_fp': to('data/simu/umi_lens/spacer/'),
                    'fastq_fn': 'umi_len_' + str(umi_len),
                },
            }
            p = simugcondi(simu_params)
            print(p.ondemand())
        return


if __name__ == "__main__":
    p = generalCondi()

    # print(p.pcrNums())

    # print(p.pcrErrs())

    # print(p.seqErrs())

    print(p.umiLens())

    # print(p.amplRates())