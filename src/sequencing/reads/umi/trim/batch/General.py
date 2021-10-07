__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from Path import to
from src.sequencing.reads.umi.trim.Template import template as umitrim
from src.sequencing.reads.simulate.dispatcher.batch.General import general as generalstarter


class general(generalstarter):

    def __init__(self, ):
        super(general, self).__init__()

    def pcrNums(self, ):
        for id, i_pcr_num in enumerate(self.pcr_nums):
            trim_params = {
                'umi_1': {
                    'len': self.umi_unit_len_fixed,
                },
                'seq_struct': 'umi_1+seq_1',
                'seq_1': {
                    'len': 100 - self.umi_unit_len_fixed,
                },
                'fastq': {
                    'path': to('data/simu/pcr_nums/'),
                    'name': 'pcr_' + str(i_pcr_num),
                    'trimmed_path': to('data/simu/pcr_nums/trimmed/'),
                    'trimmed_name': 'pcr_' + str(i_pcr_num),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)

    def pcrErrs(self, ):
        for id, i_pcr_err in enumerate(self.pcr_errs):
            trim_params = {
                'seq_struct': 'umi_1+seq_1',
                'umi_1': {
                    'len': self.umi_unit_len_fixed,
                },
                'seq_1': {
                    'len': 100 - self.umi_unit_len_fixed,
                },
                'fastq': {
                    'path': to('data/simu/pcr_errs/'),
                    'name': 'pcr_err_' + str(id),
                    'trimmed_path': to('data/simu/pcr_errs/trimmed/'),
                    'trimmed_name': 'pcr_err_' + str(id),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)

    def seqErrs(self, ):
        for id, i_seq_err in enumerate(self.seq_errs):
            trim_params = {
                'seq_struct': 'umi_1+seq_1',
                'umi_1': {
                    'len': self.umi_unit_len_fixed,
                },
                'seq_1': {
                    'len': 100 - self.umi_unit_len_fixed,
                },
                'fastq': {
                    'path': to('data/simu/seq_errs/'),
                    'name': 'seq_err_' + str(id),
                    'trimmed_path': to('data/simu/seq_errs/trimmed/'),
                    'trimmed_name': 'seq_err_' + str(id),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)

    def amplRates(self, ):
        for id, i_ampl_rate in enumerate(self.ampl_rates):
            trim_params = {
                'seq_struct': 'umi_1+seq_1',
                'umi_1': {
                    'len': self.umi_unit_len_fixed,
                },
                'seq_1': {
                    'len': 100 - self.umi_unit_len_fixed,
                },
                'fastq': {
                    'path': to('data/simu/ampl_rates/'),
                    'name': 'ampl_rate_' + str(id),
                    'trimmed_path': to('data/simu/ampl_rates/trimmed/'),
                    'trimmed_name': 'ampl_rate_' + str(id),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)

    def umiLens(self, ):
        for id, i_umi_len in enumerate(self.umi_unit_lens):
            trim_params = {
                'seq_struct': 'umi_1+seq_1',
                'umi_1': {
                    'len': i_umi_len,
                },
                'seq_1': {
                    'len': 100 - i_umi_len,
                },
                'fastq': {
                    'path': to('data/simu/umi_lens/'),
                    'name': 'umi_len_' + str(i_umi_len),
                    'trimmed_path': to('data/simu/umi_lens/trimmed/'),
                    'trimmed_name': 'umi_len_' + str(i_umi_len),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)


if __name__ == "__main__":
    p = general()

    # print(p.pcrNums())

    # print(p.pcrErrs())

    # print(p.seqErrs())

    # print(p.amplRates())

    print(p.umiLens())