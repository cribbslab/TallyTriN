__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from simreadflow.simulate.dispatcher.single.CDNA import cdna as simucdna
from Path import to


class general(object):

    def __init__(self, ):
        self.umi_unit_len_fixed = 12
        self.umi_unit_pattern_fixed = 3
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

        self.offset = '../' * 9

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

    def pcrNums(self, ):
        for id, i_pcr_num in enumerate(self.pcr_nums):
            simu_params = {
                'init_seq_setting': {
                    'cdna_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern_fixed,
                    'umi_unit_len': self.umi_unit_len_fixed,
                    'is_seed': True,
                    'primer_len': 5,
                    'is_sv_umi_lib': True if id == 0 else False,
                    'is_sv_seq_lib': True if id == 0 else False,
                    'is_sv_primer_lib': True if id == 0 else False,
                    'working_dir': to('data/cdna/trimer/simu/pcr_nums/'),
                    'umi_lib_fpn': to('data/cdna/trimer/simu/pcr_nums/umi.txt'),
                    'seq_lib_fpn': to('data/cdna/trimer/simu/pcr_nums/seq.txt'),
                    'primer_lib_fpn': to('data/cdna/trimer/simu/pcr_nums/primer.txt'),
                    'cand_pool_fpn': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
                    'cdna_fp': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
                },
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': i_pcr_num,
                'pcr_error': self.pcr_err_fixed,
                'seq_error': self.seq_err_fixed,
                'use_seed': False,
                'seed': None,
                'err_num_met': 'nbionom',
                'write': {
                    'fastq_fp': to('data/cdna/trimer/simu/pcr_nums/'),
                    'fastq_fn': 'pcr_' + str(i_pcr_num),
                }
            }
            p = simucdna(simu_params)
            print(p.ondemand())
        return

    def pcrErrs(self, ):
        for id, i_pcr_err in enumerate(self.pcr_errs):
            simu_params = {
                'init_seq_setting': {
                    'cdna_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern_fixed,
                    'umi_unit_len': self.umi_unit_len_fixed,
                    'is_seed': True,
                    'primer_len': 5,
                    'is_sv_umi_lib': True if id == 0 else False,
                    'is_sv_seq_lib': True if id == 0 else False,
                    'is_sv_primer_lib': True if id == 0 else False,
                    'working_dir': to('data/cdna/trimer/simu/pcr_errs/'),
                    'umi_lib_fpn': to('data/cdna/trimer/simu/pcr_errs/umi.txt'),
                    'seq_lib_fpn': to('data/cdna/trimer/simu/pcr_errs/seq.txt'),
                    'primer_lib_fpn': to('data/cdna/trimer/simu/pcr_errs/primer.txt'),
                    'cand_pool_fpn': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
                    'cdna_fp': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
                },
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': self.pcr_num_fixed,
                'pcr_error': i_pcr_err,
                'seq_error': self.seq_err_fixed,
                'use_seed': False,
                'seed': None,
                'err_num_met': 'nbionom',
                'write': {
                    'fastq_fp': to('data/cdna/trimer/simu/pcr_errs/'),
                    'fastq_fn': 'pcr_err_' + str(id),
                }
            }
            p = simucdna(simu_params)
            print(p.ondemand())
        return

    def seqErrs(self, ):
        for id, i_seq_err in enumerate(self.seq_errs):
            simu_params = {
                'init_seq_setting': {
                    'cdna_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern_fixed,
                    'umi_unit_len': self.umi_unit_len_fixed,
                    'is_seed': True,
                    'primer_len': 5,
                    'is_sv_umi_lib': True if id == 0 else False,
                    'is_sv_seq_lib': True if id == 0 else False,
                    'is_sv_primer_lib': True if id == 0 else False,
                    'working_dir': to('data/cdna/trimer/simu/seq_errs/'),
                    'umi_lib_fpn': to('data/cdna/trimer/simu/seq_errs/umi.txt'),
                    'seq_lib_fpn': to('data/cdna/trimer/simu/seq_errs/seq.txt'),
                    'primer_lib_fpn': to('data/cdna/trimer/simu/seq_errs/primer.txt'),
                    'cand_pool_fpn': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
                    'cdna_fp': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
                },
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': self.pcr_num_fixed,
                'pcr_error': self.pcr_err_fixed,
                'seq_error': i_seq_err,
                'use_seed': False,
                'seed': None,
                'err_num_met': 'nbionom',
                'write': {
                    'fastq_fp': to('data/cdna/trimer/simu/seq_errs/'),
                    'fastq_fn': 'seq_err_' + str(id),
                }
            }
            p = simucdna(simu_params)
            print(p.ondemand())
        return

    def umiLens(self, ):
        for id, umi_len in enumerate(self.umi_unit_lens):
            simu_params = {
                'init_seq_setting': {
                    'cdna_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern_fixed,
                    'umi_unit_len': umi_len,
                    'is_seed': True,
                    'primer_len': 5,
                    'is_sv_umi_lib': True,
                    'is_sv_seq_lib': True,
                    'is_sv_primer_lib': True,
                    'working_dir': to('data/cdna/trimer/simu/umi_lens/'),
                    'umi_lib_fpn': to('data/cdna/trimer/simu/umi_lens/umi_') + str(umi_len) + '.txt',
                    'seq_lib_fpn': to('data/cdna/trimer/simu/umi_lens/seq_') + str(umi_len) + '.txt',
                    'primer_lib_fpn': to('data/cdna/trimer/simu/umi_lens/primer_') + str(umi_len) + '.txt',
                    'cand_pool_fpn': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
                    'cdna_fp': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
                },
                'ampl_rate': self.ampl_rate_fixed,
                'pcr_num': self.pcr_num_fixed,
                'pcr_error': self.pcr_err_fixed,
                'seq_error': self.seq_err_fixed,
                'use_seed': False,
                'seed': None,
                'err_num_met': 'nbionom',
                'write': {
                    'fastq_fp': to('data/cdna/trimer/simu/umi_lens/'),
                    'fastq_fn': 'umi_len_' + str(umi_len),
                }
            }
            p = simucdna(simu_params)
            print(p.ondemand())
        return

    def amplRates(self, ):
        for id, ampl_rate in enumerate(self.ampl_rates):
            simu_params = {
                'init_seq_setting': {
                    'cdna_num': self.umi_num_fixed,
                    'umi_unit_pattern': self.umi_unit_pattern_fixed,
                    'umi_unit_len': self.umi_unit_len_fixed,
                    'is_seed': True,
                    'primer_len': 5,
                    'is_sv_umi_lib': True if id == 0 else False,
                    'is_sv_seq_lib': True if id == 0 else False,
                    'is_sv_primer_lib': True if id == 0 else False,
                    'working_dir': to('data/cdna/trimer/simu/ampl_rates/'),
                    'umi_lib_fpn': to('data/cdna/trimer/simu/ampl_rates/umi.txt'),
                    'seq_lib_fpn': to('data/cdna/trimer/simu/ampl_rates/seq.txt'),
                    'primer_lib_fpn': to('data/cdna/trimer/simu/ampl_rates/primer.txt'),
                    'cand_pool_fpn': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
                    'cdna_fp': self.offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
                },
                'ampl_rate': ampl_rate,
                'pcr_num': self.pcr_num_fixed,
                'pcr_error': self.pcr_err_fixed,
                'seq_error': self.seq_err_fixed,
                'use_seed': False,
                'seed': None,
                'err_num_met': 'nbionom',
                'write': {
                    'fastq_fp': to('data/cdna/trimer/simu/ampl_rates/'),
                    'fastq_fn': 'ampl_rate_' + str(id),
                }
            }
            p = simucdna(simu_params)
            print(p.ondemand())
        return


if __name__ == "__main__":
    p = general()

    print(p.pcrNums())

    # print(p.pcrErrs())

    # print(p.seqErrs())

    # print(p.umiLens())

    # print(p.amplRates())