__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np


class config(object):

    def __init__(self, ):
        self.permutation_num = 2

        self.umi_unit_pattern = 2
        self.umi_unit_len_fixed = 12
        # self.seq_len_fixed = 100
        self.umi_num_fixed = 50
        self.pcr_num_fixed = 16
        self.pcr_err_fixed = 1e-4
        self.seq_err_fixed = 1e-2
        self.ampl_rate_fixed = 0.80
        self.sim_thres_fixed = 3
        self.seq_sub_spl_rate = 0.3333

        self.ampl_rates = np.linspace(0.1, 1, 10)
        self.umi_unit_lens = np.arange(7, 36 + 1, 1)
        self.umi_nums = np.arange(20, 140 + 20, 20)
        self.pcr_nums = np.arange(1, 18 + 1, 1)
        # self.pcr_nums = np.arange(1, 20 + 1, 1)
        self.pcr_errs, self.seq_errs = self.errors()
        self.pcr_errs = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01]

        # self.seq_fix_errs = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2]
        # self.seq_errs = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05]
        print(self.pcr_errs, self.seq_errs)

        self.metric_vals = {
            'pcr_nums': self.pcr_nums,
            'pcr_errs': self.pcr_errs,
            # 'seq_errs': [1e-05, 2.5e-05],
            'seq_errs': self.seq_errs,
            'ampl_rates': self.ampl_rates,
            'umi_lens': self.umi_unit_lens,
        }
        self.fn_pref = {
            'pcr_nums': 'pcr_num_',
            'pcr_errs': 'pcr_err_',
            'seq_errs': 'seq_err_',
            'ampl_rates': 'ampl_rate_',
            'umi_lens': 'umi_len_',
        }

    # def errors(self, ):
    #     pcr_errs = []
    #     seq_errs = []
    #     e = 1e-5
    #     while e < 3e-1:
    #         pcr_errs.append(e)
    #         seq_errs.append(e)
    #         if 5 * e < 3e-1:
    #             pcr_errs.append(2.5 * e)
    #             pcr_errs.append(5 * e)
    #             pcr_errs.append(7.5 * e)
    #             seq_errs.append(2.5 * e)
    #             seq_errs.append(5 * e)
    #
    #             seq_errs.append(7.5 * e)
    #         e = 10 * e
    #     pcr_errs.append(0.125)
    #     seq_errs.append(0.125)
    #     pcr_errs.append(0.15)
    #     seq_errs.append(0.15)
    #     pcr_errs.append(0.2)
    #     seq_errs.append(0.2)
    #     pcr_errs.append(0.225)
    #     seq_errs.append(0.225)
    #     pcr_errs.append(0.25)
    #     seq_errs.append(0.25)
    #     pcr_errs.append(0.3)
    #     seq_errs.append(0.3)
    #     # print(pcr_errs)
    #     # print(seq_errs)
    #     return pcr_errs, seq_errs

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