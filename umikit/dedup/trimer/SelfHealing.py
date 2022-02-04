__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import textwrap
import pandas as pd
from collections import Counter
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.random.Number import number as rannum
from umikit.trim.Template import template as umitrim
from umikit.trim.Reader import reader as trimreader
from simreadflow.util.file.read.Reader import reader as gfreader
from simreadflow.util.file.write.Writer import writer as fwriter
from simreadflow.simulate.dispatcher.batch.UMI import umi as generalstarter
from simreadflow.read.similarity.distance.Hamming import hamming
from Path import to


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
            df_fpn=to('data/simu/umi/seq_errs/trimer/umi.txt')
        )[0].values
        print({i: e for i, e in enumerate(self.umi_raw)})
        self.umi_raw = pd.DataFrame.from_dict({i: e for i, e in enumerate(self.umi_raw)}, orient='index')
        print(self.umi_raw)
        self.umi_raw_monomer = self.umi_raw[0].apply(lambda x: ''.join([i[0] for i in textwrap.wrap(x, 3)])).to_dict()
        print(self.umi_raw_monomer)

    def rea(self, ):
        df_stat = pd.DataFrame()
        for id, i_seq_err in enumerate(self.seq_errs):
            read_stime = time.time()
            names, seqs, _, _ = self.rfastq().fromgz(
                fastq_path=to('data/simu/umi/seq_errs/trimer/trimmed/'),
                fastq_name='seq_err_' + str(id),
            )
            print('===>file read time: {:.3f}s'.format(time.time() - read_stime))
            df_fastq = self.trimreader.todf(names=names, seqs=seqs)
            mono_corr_stime = time.time()
            df_fastq['umi_mono_corr'] = df_fastq['umi'].apply(lambda x: self.correct(x))
            print('===>mono_corr time: {:.3f}s'.format(time.time() - mono_corr_stime))
            hm_stime = time.time()
            df_stat['umi_hm' + str(id)] = df_fastq.apply(lambda x: hamming().general(x['umi_mono_corr'], self.umi_raw_monomer[x['umi#']]), axis=1)
            print('===>hamming time: {:.3f}s'.format(time.time() - hm_stime))
            print(df_stat['umi_hm' + str(id)])
        self.fwriter.generic(df=df_stat, sv_fpn=to('data/simu/umi/seq_errs/trimer/trimmed/dasd.txt'), df_sep='\t')
        return

    def correct(self, umi):
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        # umi_trimers = textwrap.wrap(umi, 3)
        t = []
        for umi_trimer in umi_trimers:
            s = set(umi_trimer)
            if len(s) == 3:
                rand_index = self.rannum.uniform(low=0, high=3, num=1, use_seed=False)[0]
                t.append(umi_trimer[rand_index])
            elif len(s) == 2:
                sdict = {umi_trimer.count(i): i for i in s}
                t.append(sdict[2])
            else:
                t.append(umi_trimer[0])
        return ''.join(t)

    def correct_deprecate(self, umi):
        vernier = [i for i in range(36) if i % 3 == 0]
        umi_trimers = [umi[v: v+3] for v in vernier]
        # umi_trimers = textwrap.wrap(umi, 3)
        t = []
        for umi_trimer in umi_trimers:
            s = Counter(umi_trimer).most_common()
            if len(s) == 3:
                rand_index = self.rannum.uniform(low=0, high=3, num=1, use_seed=False)[0]
                t.append(s[rand_index][0])
            elif len(s) == 2:
                t.append(s[0][0])
            else:
                t.append(umi_trimer[0])
        return ''.join(t)

    def stat(self, x):
        return hamming().general(x['umi_mono_corr'], self.umi_raw_monomer[x['umi#']])


if __name__ == "__main__":
    p = selfHealing()

    print(p.rea())