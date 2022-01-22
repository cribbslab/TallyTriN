__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time


class reader(object):

    def __init__(self, ):
        pass

    def todf(self, seqs, names):
        import pandas as pd
        umi_df_stime = time.time()
        df_fastq = pd.DataFrame(seqs, columns=['seq_raw'])
        df_fastq['name'] = names
        df_fastq['umi'] = df_fastq['name'].apply(lambda x: x.split('_')[1])
        df_fastq['umi#'] = df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[0])
        df_fastq['umi_src'] = df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[1])
        # df_fastq['umi_pcr#'] = df_fastq['name'].apply(lambda x: self.pcrnum(x))
        df_fastq['umi#'] = df_fastq['umi#'].astype(int)
        # print(df_fastq)
        # print(df_fastq['name'])
        # print(df_fastq['umi'])
        # print(df_fastq['umi#'])
        # print(df_fastq['umi_src'])
        print('===>Trimmed UMIs to df time: {:.3f}s'.format(time.time() - umi_df_stime))
        return df_fastq

    def pcrnum(self, x):
        c = x.split('_')[0].split('-')
        if c[1] == 'init':
            return -1
        else:
            return c[2]