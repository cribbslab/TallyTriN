__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import textwrap
import numpy as np
import pandas as pd
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.sequence.fastq.Write import write as wfastq
from simreadflow.util.random.Number import number as rannum
from simreadflow.util.file.read.Reader import reader as gfreader
from simreadflow.util.file.write.Writer import writer as fwriter
from umikit.trim.Template import template as umitrim
from umikit.trim.Reader import reader as trimreader
from umikit.dedup.dimer.pipeline import Config
from Path import to


class selfHealing(Config.config):

    def __init__(self, fastq_fp, cat):
        super(selfHealing, self).__init__()
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
                    fastq_path=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/trimmed/',
                    fastq_name=self.cat + '_' + str(id),
                    method='pyfastx',
                )
                print('permutation {}, No.{} with criterion: {}'.format(i_pn, id, i_seq_err))
                # df = self.trimreader.todf(names=names, seqs=seqs)
                df = self.trimreader.todfFromTree(names=names, seqs=seqs)
                # print(df)
                df['origin_info'] = df.apply(lambda x: str(x['umi#']) + '_' + x['umi_src'], axis=1)
                mono_corr_stime = time.time()
                df['umi_monos'] = df['umi'].apply(lambda x: self.split(x, num=20))
                df['umi_correct'] = df['umi'].apply(lambda x: self.correct(x))
                df['umi_mark'] = df['umi_monos'].apply(lambda x: self.marker(x))
                df['umi_l'] = df['umi_monos'].apply(lambda x: x.split(';')[0])
                df['umi_r'] = df['umi_monos'].apply(lambda x: x.split(';')[1])
                df['to_fas'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_l'], axis=1)
                df['to_fas_correct'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_correct'], axis=1)
                df_fastq_cp = df.copy()
                df_fastq_cp['to_fas'] = df_fastq_cp.apply(lambda x: x['origin_info'] + '_' + x['umi_r'], axis=1)
                df_umi_extra = df.loc[df['umi_mark'] == 1]
                df_umi_extra['to_fas'] = df_umi_extra.apply(lambda x: x['origin_info'] + '_' + x['umi_r'], axis=1)
                # print(df_umi_extra['to_fas'])
                df_merge = pd.concat(
                    [df[['seq_raw', 'to_fas']], df_fastq_cp[['seq_raw', 'to_fas']]],
                    axis=0,
                ).reset_index(drop=True)
                print(df_merge)
                df_extra = pd.concat(
                    [df[['seq_raw', 'to_fas']], df_umi_extra[['seq_raw', 'to_fas']]],
                    axis=0,
                ).reset_index(drop=True)
                print(df_extra)
                self.wfastq.togz(
                    list_2d=df[['seq_raw', 'to_fas']].values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/ref/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df[['seq_raw', 'to_fas_correct']].values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/correct/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df_merge.values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/double/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df_extra.values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/bipartite/',
                    fn=self.cat + '_' + str(id),
                )
                print('===>monos time: {:.3f}s'.format(time.time() - mono_corr_stime))
        return

    def split(self, umi, num=24):
        vernier = [i for i in range(num) if i % 2 == 0]
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

    def marker(self, umi):
        tt = umi.split(';')
        if tt[0] != tt[1]:
            return 1
        else:
            return 0


if __name__ == "__main__":
    p = selfHealing(
        fastq_fp=to('data/simu/dimer/tree250/'),
        cat='seq_err',
    )
    print(p.rea())