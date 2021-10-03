__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from Path import to
from src.sequencing.reads.umi.trim.Template import template as umitrim
from src.sequencing.reads.simulate.dispatcher.batch.UMI import umi as generalstarter


class general(generalstarter):

    def __init__(self, ):
        super(general, self).__init__()

    def pcrNums(self, ):
        pass

    def pcrErrs(self, ):
        pass

    def seqErrs(self, ):
        for id, i_seq_err in enumerate(self.seq_errs):
            trim_params = {
                'seq_struct': 'umi_1',
                'umi_1': {
                    'len': self.umi_unit_len_fixed * self.umi_unit_pattern,
                },
                'fastq': {
                    'path': to('data/simu/umi/seq_errs/trimer/'),
                    'name': 'seq_err_' + str(id),
                    'trimmed_path': to('data/simu/umi/seq_errs/trimer/trimmed/'),
                    'trimmed_name': 'seq_err_' + str(id),
                },
            }
            umitrim_parser = umitrim(trim_params)
            df = umitrim_parser.todf()
            umitrim_parser.togz(df)

    def amplRates(self, ):
        pass

    def umiLens(self, ):
        pass


if __name__ == "__main__":
    p = general()

    # print(p.pcrNums())

    # print(p.pcrErrs())

    print(p.seqErrs())

    # print(p.amplRates())

    # print(p.umiLens())