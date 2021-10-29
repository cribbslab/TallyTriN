__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
import time
from src.util.sequence.fastq.Read import read as rfastq
from src.util.random.Number import number as rannum
from src.sequencing.reads.umi.trim.Template import template as umitrim
from src.sequencing.reads.umi.trim.Reader import reader as trimreader
from Path import to
from src.util.file.read.Reader import reader as gfreader
from src.sequencing.reads.simulate.dispatcher.batch.UMI import umi as generalstarter
import textwrap
from src.sequencing.reads.similarity.distance.Hamming import hamming
from src.util.file.write.Writer import writer as fwriter


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
            df_fpn=to('data/simu/umi/seq_errs/dimer/umi.txt')
        )[0].values
        print({i: e for i, e in enumerate(self.umi_raw)})
        self.umi_raw = pd.DataFrame.from_dict({i: e for i, e in enumerate(self.umi_raw)}, orient='index')
        self.umi_mono_raw = self.umi_raw[0].apply(lambda x: ''.join([i[0] for i in textwrap.wrap(x, 2)])).to_dict()
        print(self.umi_mono_raw)

    def rea(self, ):
        df_stat = pd.DataFrame()
        for id, i_seq_err in enumerate(self.seq_errs):
            read_stime = time.time()
            names, seqs, _, _ = self.rfastq().fromgz(
                fastq_path=to('data/simu/umi/seq_errs/dimer/trimmed/'),
                fastq_name='seq_err_' + str(id),
            )
            print('===>file read time: {:.3f}s'.format(time.time() - read_stime))
            df_fastq = self.trimreader.todf(names=names, seqs=seqs)
            mono_corr_stime = time.time()
            df_fastq['umi_mono_corr'] = df_fastq['umi'].apply(lambda x: self.correct(x))
            print('===>mono_corr time: {:.3f}s'.format(time.time() - mono_corr_stime))
            hm_stime = time.time()
            df_stat['umi_hm' + str(id)] = df_fastq.apply(lambda x: hamming().general(x['umi_mono_corr'], self.umi_mono_raw[x['umi#']]), axis=1)
            print('===>hamming time: {:.3f}s'.format(time.time() - hm_stime))
            tui = df_stat['umi_hm' + str(id)].values
            print(len(tui[tui != 0]))
        self.fwriter.generic(df=df_stat, sv_fpn=to('data/simu/umi/seq_errs/dimer/trimmed/dasd.txt'), df_sep='\t')
        return

    def correct(self, umi):
        umi_dimers = textwrap.wrap(umi, 2)
        t = []
        # print(umi_dimers)
        for umi_dimer in umi_dimers:
            if umi_dimer[0] != umi_dimer[1]:
                rand_index = np.random.randint(low=0, high=2, size=1)[0]
                # print(rand_index)
                t.append(umi_dimer[rand_index])
            else:
                t.append(umi_dimer[0])
        return ''.join(t)


if __name__ == "__main__":
    p = selfHealing()
    print(p.rea())