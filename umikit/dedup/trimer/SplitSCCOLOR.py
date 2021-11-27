__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import textwrap
import pandas as pd
from Path import to
from collections import Counter
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.sequence.fastq.Write import write as wfastq
from simreadflow.util.random.Number import number as rannum
from umikit.trim.Template import template as umitrim
from umikit.trim.Reader import reader as trimreader
from simreadflow.util.file.read.Reader import reader as gfreader
from simreadflow.util.file.write.Writer import writer as fwriter
from umikit.dedup.trimer.pipeline import Config
from simreadflow.read.similarity.distance.Hamming import hamming


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
        # for i_pn in [2]:
            for id, i_seq_err in enumerate(self.seq_errs):
                read_stime = time.time()
                print(self.fastq_fp  + '/permute_' + str(i_pn) + '/trimmed/' + self.cat + '_' + str(id))
                names, seqs, _, _ = self.rfastq.fromgz(
                    fastq_path=self.fastq_fp  + '/permute_' + str(i_pn) + '/trimmed/',
                    fastq_name=self.cat + '_' + str(id),
                    method='pyfastx',
                )
                print('===>file read time: {:.3f}s'.format(time.time() - read_stime))
                df_fastq = self.trimreader.todf(names=names, seqs=seqs)
                df_fastq['origin_info'] = df_fastq['name'].apply(lambda x: x.split('_')[0])
                mono_corr_stime = time.time()

                # df_fastq['umi_monos'] = df_fastqNo3['umi'].apply(lambda x: self.split(x))
                df_fastq['umi_monos'] = df_fastq['umi'].apply(lambda x: self.split(x))
                print(df_fastq['umi_monos'])
                df_fastq['umi_raw'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[0])
                df_fastq['umi_l'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[1])
                df_fastq['umi_m'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[2])
                df_fastq['umi_r'] = df_fastq['umi_monos'].apply(lambda x: x.split(';')[3])

                df_fastq['to_fas'] = df_fastq.apply(lambda x: x['origin_info'] + '_' + x['umi_raw'], axis=1)

                df_fastq['umi_mark'] = df_fastq['umi'].apply(lambda x: self.marker(x))
                df_fastqYes3 = df_fastq.loc[df_fastq['umi_mark'] == 1]
                df_fastqYes3_cp = df_fastqYes3.copy()
                df_fastqYes3 = df_fastqYes3.drop('umi_r', 1)
                df_fastqYes3_cp = df_fastqYes3_cp.drop('umi_m', 1)
                df_fastqYes3 = df_fastqYes3.rename(columns={"umi_m": "umi_bi"})
                df_fastqYes3_cp = df_fastqYes3_cp.rename(columns={"umi_r": "umi_bi"})
                df_Yes3 = pd.concat([df_fastqYes3, df_fastqYes3_cp], axis=0).reset_index(drop=True)
                # print(df_fastqYes3['umi_bi'])
                # print(df_fastqYes3_cp['umi_bi'])
                # print(df_Yes3['umi_bi'])
                print('jianfeng')
                df_fastqNo3 = df_fastq.loc[df_fastq['umi_mark'] != 1]
                df_fastqNo3['to_fas'] = df_fastqNo3.apply(lambda x: x['origin_info'] + '_' + x['umi_l'], axis=1)
                if df_Yes3.empty:
                    df = df_fastqNo3[['seq_raw', 'to_fas']].reset_index(drop=True)
                else:
                    df_Yes3['to_fas'] = df_Yes3.apply(lambda x: x['origin_info'] + '_' + x['umi_bi'], axis=1)
                    df = pd.concat([df_Yes3[['seq_raw', 'to_fas']], df_fastqNo3[['seq_raw', 'to_fas']]], axis=0).reset_index(drop=True)
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

                # print(df_fastq['umi_monosa'])
                print('===>mono_corr time: {:.3f}s'.format(time.time() - mono_corr_stime))
        return

    def split(self, umi):
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        ref = []
        l = []
        m = []
        r = []
        for umi_trimer in umi_trimers:
            ref.append(umi_trimer[0])
            s = set(umi_trimer)
            if len(s) == 3:
                l.append(umi_trimer[0])
                m.append(umi_trimer[1])
                r.append(umi_trimer[2])
            elif len(s) == 2:
                sdict = {umi_trimer.count(i): i for i in s}
                l.append(sdict[2])
                m.append(sdict[2])
                r.append(sdict[2])
            else:
                l.append(umi_trimer[0])
                m.append(umi_trimer[0])
                r.append(umi_trimer[0])
        ref = ''.join(ref)
        l = ''.join(l)
        m = ''.join(m)
        r = ''.join(r)
        return ref + ';' + l + ';' + m + ';' + r

    def marker(self, umi):
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        slens = [len(set(umi_trimer)) for umi_trimer in umi_trimers]
        if 3 in slens:
            return 1
        else:
            return 0


if __name__ == "__main__":
    p = selfHealing(
        fastq_fp=to('data/simu/trimer/pcr8/seq_errs/'),
        cat='seq_err',
    )

    print(p.rea())