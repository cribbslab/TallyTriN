__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import pandas as pd
import time
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.random.Number import number as rannum
from umikit.trim.Template import template as umitrim
from umikit.trim.Reader import reader as trimreader
from Path import to
from simreadflow.util.file.read.Reader import reader as gfreader
from simreadflow.simulate.dispatcher.batch.UMI import umi as generalstarter
from simreadflow.read.similarity.distance.Hamming import hamming
from simreadflow.util.file.write.Writer import writer as fwriter


class selfHealing(generalstarter):

    def __init__(self, ):
        super(selfHealing, self).__init__()
        self.umitrim = umitrim
        self.gfreader = gfreader()
        self.rfastq = rfastq
        self.trimreader = trimreader()
        self.rannum = rannum()
        self.fwriter = fwriter()
        self.umi_raw = self.gfreader.generic(
            df_fpn=to('data/simu/umi/seq_errs/monomer/umi.txt')
        )[0].values
        self.umi_raw_monomer = {i: e for i, e in enumerate(self.umi_raw)}
        print(self.umi_raw_monomer)

    def rea(self, ):
        df_stat = pd.DataFrame()
        for id, i_seq_err in enumerate(self.seq_errs):
            read_stime = time.time()
            names, seqs, _, _ = self.rfastq().fromgz(
                fastq_path=to('data/simu/umi/seq_errs/monomer/trimmed/'),
                fastq_name='seq_err_' + str(id),
            )
            print('===>file read time: {:.3f}s'.format(time.time() - read_stime))
            df_fastq = self.trimreader.todf(names=names, seqs=seqs)
            hm_stime = time.time()
            df_stat['umi_hm' + str(id)] = df_fastq.apply(lambda x: hamming().general(x['umi'], self.umi_raw_monomer[x['umi#']]), axis=1)
            print('===>hamming time: {:.3f}s'.format(time.time() - hm_stime))
            print(df_stat['umi_hm' + str(id)])
        self.fwriter.generic(df=df_stat, sv_fpn=to('data/simu/umi/seq_errs/monomer/trimmed/dasd.txt'), df_sep='\t')
        return

    def correct(self, umi):
        pass

    def stat(self, x):
        return hamming().general(x['umi_mono_corr'], self.umi_raw_monomer[x['umi#']])


if __name__ == "__main__":
    p = selfHealing()

    print(p.rea())