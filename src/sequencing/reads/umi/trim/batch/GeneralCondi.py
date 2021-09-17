__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import numpy as np
from Path import to
from src.sequencing.reads.umi.trim.Template import template as umitrim
from src.sequencing.reads.simulate.starter.batch.General import general as generalstarter


class general(generalstarter):

    def __init__(self, ):
        super(general, self).__init__()

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
                    'path': to('data/simu/umi_lens/spacer/'),
                    'name': 'umi_len_' + str(i_umi_len),
                    'trimmed_path': to('data/simu/umi_lens/spacer/trimmed1/'),
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