__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import pandas as pd
import time
from src.util.sequence.fastq.Read import read as rfastq
from src.sequencing.reads.umi.trim.Template import template as umitrim
from Path import to
from src.util.file.read.Reader import reader as gfreader
from collections import Counter
from src.sequencing.reads.simulate.dispatcher.batch.UMI import umi as generalstarter


class selfHealing(generalstarter):

    def __init__(self, ):
        super(selfHealing, self).__init__()
        self.umitrim = umitrim
        self.gfreader = gfreader()
        self.rfastq = rfastq

    def rea(self, ):
        umi = self.gfreader.generic(
            df_fpn=to('data/simu/umi/seq_errs/monomer/umi.txt')
        )[0].values
        umi = {i: e for i, e in enumerate(umi)}
        print(umi)
        for id, i_seq_err in enumerate(self.seq_errs):
            read_stime = time.time()
            names, seqs, _, _ = self.rfastq().fromgz(
                fastq_path=to('data/simu/umi/seq_errs/monomer/trimmed/'),
                fastq_name='seq_err_' + str(id),
            )

            print('--->file read time: {:.3f}s'.format(time.time() - read_stime))
            umi_df_stime = time.time()
            df_fastq = pd.DataFrame(seqs, columns=['seq_raw'])
            df_fastq['name'] = names
            df_fastq['umi'] = df_fastq['name'].apply(lambda x: x.split('_')[1])
            df_fastq['umi#'] = df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[0])
            df_fastq['umi_src'] = df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[1])
            # df_fastq['umi_pcr#'] = df_fastq['name'].apply(lambda x: pcrnum(x))
            print(df_fastq)
            print(df_fastq['name'])
            print(df_fastq['umi'])
            print(df_fastq['umi#'])
            print(df_fastq['umi_src'])
            print('--->umi to df time: {:.3f}s'.format(time.time() - umi_df_stime))

        return


if __name__ == "__main__":
    p = selfHealing()
    print(p.rea())