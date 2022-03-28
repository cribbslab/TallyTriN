__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import sys
import numpy as np
from simreadflow.simulate.dispatcher.single.UMIDouble import umiDouble as simugeneral
from Path import to


class umiDouble(object):

    def __init__(self, working_dir):
        # ### /*** block. general ***/
        self.permutation_num = 1000

        self.working_dir = working_dir
        self.umi_unit_pattern = 3
        self.umi_unit_len_fixed = 10
        self.umi_num_fixed = 100
        self.pcr_num_fixed = 8
        self.pcr_err_fixed = 1e-3
        self.seq_err_fixed = 1e-3
        self.ampl_rate_fixed = 0.85
        self.sim_thres_fixed = 3
        self.seq_sub_spl_rate = 1

        self.ampl_rates = np.linspace(0.1, 1, 10)
        self.umi_unit_lens = np.arange(6, 36 + 1, 1)
        self.umi_nums = np.arange(20, 140 + 20, 20)
        self.pcr_nums = np.arange(1, 20 + 1, 1)
        self.pcr_errs, self.seq_errs = self.errors()
        print(self.pcr_errs)
        print(self.seq_errs)

        # ### /*** block. bulk ***/
        # self.permutation_num = 10

        # self.umi_unit_pattern = 3
        # self.umi_unit_len_fixed = 12
        # self.umi_num_fixed = 100
        # self.pcr_num_fixed = 8
        # self.pcr_err_fixed = 1e-3
        # self.seq_err_fixed = 1e-3
        # self.ampl_rate_fixed = 0.85
        # self.sim_thres_fixed = 3
        # self.seq_sub_spl_rate = 1

        # self.ampl_rates = np.linspace(0.1, 1, 10)
        # self.umi_unit_lens = np.arange(6, 36 + 1, 1)
        # self.umi_nums = np.arange(20, 140 + 20, 20)
        # self.pcr_nums = np.arange(1, 20 + 1, 1)
        # self.pcr_errs, self.seq_errs = self.errors()
        # self.pcr_errs = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01]
        # print(self.pcr_errs, self.seq_errs)

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
        pcr_errs.append(0.125)
        seq_errs.append(0.125)
        pcr_errs.append(0.15)
        seq_errs.append(0.15)
        pcr_errs.append(0.2)
        seq_errs.append(0.2)
        pcr_errs.append(0.225)
        seq_errs.append(0.225)
        pcr_errs.append(0.25)
        seq_errs.append(0.25)
        pcr_errs.append(0.3)
        seq_errs.append(0.3)
        # print(pcr_errs)
        # print(seq_errs)
        return pcr_errs, seq_errs

    def seqErrs(self, ):
        sys.stdout = open(self.working_dir + 'log.txt', 'w')
        for pn in range(self.permutation_num):
            print(pn)
            simu_params = {
                'init_seq_setting': {
                    'seq_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern,
                    'umi_unit_len': self.umi_unit_len_fixed,
                    'is_seed': True,

                    # # ### /*** block. trimer ***/
                    # 'working_dir': to('data/simu/transloc/trimer/single_read/seq_errs/permute_') + str(pn) + '/',
                    # 'is_sv_umi_lib': True,
                    # 'umi_lib_fpn': to('data/simu/transloc/trimer/single_read/seq_errs/permute_') + str(pn) + '/umi.txt',

                    # ### /*** block. trimer ***/
                    'is_sv_umi_lib': True,
                    'working_dir': self.working_dir + 'seq_errs/permute_' + str(pn) + '/',
                    'umi_lib_fpn': self.working_dir + 'seq_errs/permute_' + str(pn) + '/umi.txt',

                    'condis': ['umi'],
                    'sim_thres': self.sim_thres_fixed,
                    'permutation': pn,
                },
                'transloc_rate': 0.02,
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': self.pcr_num_fixed,
                'err_num_met': 'nbinomial',
                'pcr_error': self.pcr_err_fixed,
                'seq_errors': self.seq_errs,
                'seq_sub_spl_rate': self.seq_sub_spl_rate,
                'use_seed': False,
                'seed': None,
                'write': {
                    'fastq_fp': self.working_dir + 'seq_errs/permute_' + str(pn) + '/',
                    'fastq_fn': '',
                }
            }
            p = simugeneral(simu_params)
            print(p.ondemandSeqErrs())
        sys.stdout.close()
        return


if __name__ == "__main__":
    p = umiDouble(
        working_dir=to('data/simu/transloc/trimer/single_read/1000/')
    )

    # print(p.pcrNums())

    # print(p.pcrErrs())

    print(p.seqErrs())

    # print(p.umiLens())

    # print(p.amplRates())