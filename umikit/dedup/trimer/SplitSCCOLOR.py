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
                print('permutation {}, No.{} with criterion: {}'.format(i_pn, id, i_seq_err))
                names, seqs, _, _ = self.rfastq.fromgz(
                    fastq_path=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/trimmed/',
                    fastq_name=self.cat + '_' + str(id),
                    method='pyfastx',
                )
                # df = self.trimreader.todf(names=names, seqs=seqs)
                df = self.trimreader.todfFromTree(names=names, seqs=seqs)
                # print(df)
                df['origin_info'] = df['name'].apply(lambda x: x.split('_')[0])
                mono_corr_stime = time.time()

                df['umi_extract'] = df['umi'].apply(lambda x: self.extract(x, num=30))
                df['umi_l_extract'] = df['umi_extract'].apply(lambda x: x.split(';')[0])
                df['umi_m_extract'] = df['umi_extract'].apply(lambda x: x.split(';')[1])
                df['umi_r_extract'] = df['umi_extract'].apply(lambda x: x.split(';')[2])
                df_fastq_cp_l = df.copy()
                df_fastq_cp_m = df.copy()
                df_fastq_cp_r = df.copy()
                df_fastq_cp_l['to_fas'] = df_fastq_cp_l.apply(lambda x: x['origin_info'] + '_' + x['umi_l_extract'], axis=1)
                df_fastq_cp_m['to_fas'] = df_fastq_cp_m.apply(lambda x: x['origin_info'] + '_' + x['umi_m_extract'], axis=1)
                df_fastq_cp_r['to_fas'] = df_fastq_cp_r.apply(lambda x: x['origin_info'] + '_' + x['umi_r_extract'], axis=1)
                df_lmr = pd.concat(
                    [
                        df_fastq_cp_l[['seq_raw', 'to_fas']],
                        df_fastq_cp_m[['seq_raw', 'to_fas']],
                        df_fastq_cp_r[['seq_raw', 'to_fas']],
                    ],
                    axis=0,
                ).reset_index(drop=True)

                df['umi_monos'] = df['umi'].apply(lambda x: self.split(x, num=30))
                df['umi_ref'] = df['umi_monos'].apply(lambda x: x.split(';')[0])
                df['umi_l'] = df['umi_monos'].apply(lambda x: x.split(';')[1])
                df['umi_m'] = df['umi_monos'].apply(lambda x: x.split(';')[2])
                df['umi_r'] = df['umi_monos'].apply(lambda x: x.split(';')[3])

                df['to_fas'] = df.apply(lambda x: x['origin_info'] + '_' + x['umi_ref'], axis=1)
                df['umi_mark'] = df['umi'].apply(lambda x: self.marker(x, num=30))
                # print(df['umi_mark'])
                df_3differ = df.loc[df['umi_mark'] == 1]
                df_3differ_cp = df_3differ.copy()
                df_3differ = df_3differ.drop('umi_r', 1)
                df_3differ_cp = df_3differ_cp.drop('umi_m', 1)
                df_3differ = df_3differ.rename(columns={"umi_m": "umi_bi"})
                df_3differ_cp = df_3differ_cp.rename(columns={"umi_r": "umi_bi"})

                df_umi_3differ = pd.concat([df_3differ, df_3differ_cp], axis=0).reset_index(drop=True)
                print(df_umi_3differ)

                df_umi_3differ['to_fas'] = df_umi_3differ.apply(lambda x: x['origin_info'] + '_' + x['umi_bi'], axis=1)
                print(df_3differ['umi_bi'])
                # print(df_umi_3differ['umi_bi'])
                df_not3differ = df.loc[df['umi_mark'] != 1]
                df_not3differ['to_fas'] = df_not3differ.apply(lambda x: x['origin_info'] + '_' + x['umi_l'], axis=1)
                if df_umi_3differ.empty:
                    df_merge = df_not3differ[['seq_raw', 'to_fas']].reset_index(drop=True)
                else:
                    df_merge = pd.concat(
                        [df_umi_3differ[['seq_raw', 'to_fas']], df_not3differ[['seq_raw', 'to_fas']]],
                        axis=0,
                    ).reset_index(drop=True)
                print(df_merge)
                self.wfastq.togz(
                    list_2d=df_lmr[['seq_raw', 'to_fas']].values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/lmr/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df[['seq_raw', 'to_fas']].values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/ref/',
                    fn=self.cat + '_' + str(id),
                )
                self.wfastq.togz(
                    list_2d=df_merge.values,
                    sv_fp=self.fastq_fp + 'seq_errs/permute_' + str(i_pn) + '/bipartite/',
                    fn=self.cat + '_' + str(id),
                )
                # print(df['umi_monosa'])
                print('===>getting it done with time: {:.3f}s'.format(time.time() - mono_corr_stime))
        return

    def split(self, umi, num=36):
        vernier = [i for i in range(num) if i % 3 == 0]
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

    def extract(self, umi, num=36):
        vernier = [i for i in range(num) if i % 3 == 0]
        umi_trimers = [umi[v: v + 3] for v in vernier]
        l = []
        m = []
        r = []
        for umi_trimer in umi_trimers:
            l.append(umi_trimer[0])
            m.append(umi_trimer[1])
            r.append(umi_trimer[2])
        l = ''.join(l)
        m = ''.join(m)
        r = ''.join(r)
        return l + ';' + m + ';' + r

    def marker(self, umi, num=36):
        vernier = [i for i in range(num) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        slens = [len(set(umi_trimer)) for umi_trimer in umi_trimers]
        if 3 in slens:
            return 1
        else:
            return 0


if __name__ == "__main__":
    p = selfHealing(
        fastq_fp=to('data/simu/trimer/tree1000/'),
        cat='seq_err',
    )
    print(p.rea())