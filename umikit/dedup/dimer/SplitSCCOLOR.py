__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
import time
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.sequence.fastq.Write import write as wfastq
from simreadflow.util.random.Number import number as rannum
from umikit.trim.Template import template as umitrim
from umikit.trim.Reader import reader as trimreader
from Path import to
from simreadflow.util.file.read.Reader import reader as gfreader
from umikit.dedup.dimer.pipeline import Config
import textwrap
from simreadflow.read.similarity.distance.Hamming import hamming
from simreadflow.util.file.write.Writer import writer as fwriter


class selfHealing(Config.config):

    def __init__(self, umi_ref_fp, fastq_fp, cat):
        super(selfHealing, self).__init__()
        self.umi_ref_fp = umi_ref_fp
        self.fastq_fp = fastq_fp
        self.cat = cat
        self.umitrim = umitrim
        self.gfreader = gfreader()
        self.rfastq = rfastq()
        self.wfastq = wfastq()
        self.trimreader = trimreader()
        self.rannum = rannum()
        self.fwriter = fwriter()

    def rea(self, ):
        for i_pn in range(self.permutation_num):
            for id, i_seq_err in enumerate(self.seq_errs):
                read_stime = time.time()
                names, seqs, _, _ = self.rfastq.fromgz(
                    fastq_path=self.fastq_fp  + '/permute_' + str(i_pn) + '/trimmed/',
                    fastq_name=self.cat + '_' + str(id),
                    method='pyfastx',
                )
                print('===>file read time: {:.3f}s'.format(time.time() - read_stime))
                df_fastq = self.trimreader.todf(names=names, seqs=seqs)
                df_fastq['origin_info'] = df_fastq['name'].apply(lambda x: x.split('_')[0])
                mono_corr_stime = time.time()
                df_fastq['umi_monos'] = df_fastq['umi'].apply(lambda x: self.split(x))
                df_fastq['umi_mark'] = df_fastq['umi_monos'].apply(lambda x: self.marker(x))
                df_fastq['umi_l'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[0])
                df_fastq['umi_r'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[1])
                print(df_fastq['umi_r'])
                print(df_fastq.columns)
                df_fastq['to_fas'] = df_fastq.apply(lambda x: x['origin_info'] + '_' + x['umi_l'], axis=1)
                print(df_fastq['to_fas'])
                df_umi_extra = df_fastq.loc[df_fastq['umi_mark'] == 1]
                df_umi_extra['to_fas'] = df_umi_extra.apply(lambda x: x['origin_info'] + '_' + x['umi_r'], axis=1)
                print(df_umi_extra['to_fas'])
                df = pd.concat([df_fastq[['seq_raw', 'to_fas']], df_umi_extra[['seq_raw', 'to_fas']]], axis=0).reset_index(drop=True)
                print(df)
                self.wfastq.togz(
                    list_2d=df_fastq[['seq_raw', 'to_fas']].values,
                    sv_fp=self.fastq_fp + '/permute_' + str(i_pn) + '/ref/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df.values,
                    sv_fp=self.fastq_fp + '/permute_' + str(i_pn) + '/bipartite/',
                    fn=self.cat + '_' + str(id),
                )
                print('===>monos time: {:.3f}s'.format(time.time() - mono_corr_stime))
        return

    def split(self, umi):
        vernier = [i for i in range(24) if i % 2 == 0]
        umi_dimers = [umi[v: v + 2] for v in vernier]
        l = []
        r = []
        for umi_dimer in umi_dimers:
            l.append(umi_dimer[0])
            r.append(umi_dimer[1])
        l = ''.join(l)
        r = ''.join(r)
        # print(l)
        # print(r)
        return l + ';' + r

    def marker(self, umi):
        tt = umi.split(';')
        if tt[0] != tt[1]:
            return 1
        else:
            return 0


if __name__ == "__main__":
    p = selfHealing(
        umi_ref_fp=to('data/simu/dimer/pcr8/seq_errs/'),
        fastq_fp=to('data/simu/dimer/pcr8/seq_errs/'),
        cat='seq_err',
    )
    print(p.rea())