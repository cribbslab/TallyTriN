__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.util.sequence.fastq.Read import read as rfastq
from src.util.sequence.fastq.Write import write as wfastq
from src.sequencing.reads.polya.identify.PatternMatching import patternMatching as pm
from Path import to
import time


class dashboard():

    def __init__(self, *args):
        self.args = args[0]
        self.rfastq = rfastq
        self.wfastq = wfastq
        self.pm = pm()

    def pmatching(self, ):
        names, seqs, _, _ = self.rfastq().fromgz(
            fastq_path=self.args['fastq']['path'],
            fastq_name=self.args['fastq']['name'],
            method='pyfastx',
        )
        df = pd.DataFrame(seqs, columns=['seq_raw'])
        df['name'] = names
        print(df)
        eds = self.args['ed']
        polyT_ed = {}
        for ed in eds:
            stime = time.time()
            cond = '{s<=' + str(ed) + '}'
            df['polyT_num_ed' + str(ed)] = df.apply(lambda x: self.pm.regexfindall(
                cond=cond,
                query=x['seq_raw'],
                std_pattern=self.args['std_pattern'],
            ), axis=1)
            polyT_ed[str(ed)] = df['polyT_num_ed' + str(ed)].values.tolist()
            print(df.loc[df['polyT_num_ed' + str(ed)] != 0].shape)
            print('ed' + str(ed), time.time() - stime)
        c = ['palevioletred', 'slategray', 'grey']
        ccont = 0
        for k, v in polyT_ed.items():
            sns.distplot(
                v,
                kde=True,
                hist=False,
                kde_kws={
                    "color": c[ccont],
                    "lw": 3.5,
                    "label": 'ed' + str(k),
                    "alpha": 0.75,
                },
            )
            ccont += 1
        plt.xlabel('Number of potential polyT patterns', fontsize=12)
        plt.ylabel('Density', fontsize=12)
        plt.legend(fontsize=12)
        plt.show()
        return

    def pmatchin1(self, ):
        names, seqs, _, _ = self.rfastq().fromgz(
            fastq_path=self.args['fastq']['path'],
            fastq_name=self.args['fastq']['name'],
            method='pyfastx',
        )
        df = pd.DataFrame(seqs, columns=['seq_raw'])
        df['name'] = names
        df['seq_len'] = df.apply(lambda x: len(x['seq_raw']), axis=1)
        print(df)
        # sns.distplot(
        #     df['seq_len'].values.tolist(),
        #     kde=True,
        #     hist=False,
        #     kde_kws={
        #         "color": 'teal',
        #         "lw": 1.5,
        #         "alpha": 0.75,
        #     },
        # )
        # plt.xlabel('Length of reads', fontsize=12)
        # plt.ylabel('Density', fontsize=12)
        # plt.xticks(np.linspace(0, 300000,20).astype(np.int), np.linspace(0, 300000,20).astype(np.int), rotation=45, fontsize=7)
        # plt.legend(fontsize=12)
        # plt.show()
        eds = self.args['ed']
        polyT_ed = {}
        for ed in eds:
            stime = time.time()
            cond = '{s<=' + str(ed) + '}'
            df['polyT_num_ed' + str(ed)] = df.apply(lambda x: self.pm.regexfindall(
                cond=cond,
                query=x['seq_raw'],
                std_pattern=self.args['std_pattern'],
            ), axis=1)
            t = df.loc[df['polyT_num_ed' + str(ed)] != 0]
            polyT_ed[str(ed)] = t['seq_len'].values.tolist()
            print(df.loc[df['polyT_num_ed' + str(ed)] != 0].shape)
            print('ed' + str(ed), time.time() - stime)
        c = ['palevioletred', 'slategray', 'grey']
        ccont = 0
        for k, v in polyT_ed.items():
            sns.distplot(
                v,
                kde=True,
                hist=False,
                kde_kws={
                    "color": c[ccont],
                    "lw": 2.5,
                    "label": 'edit distance ' + str(k),
                    "alpha": 0.75,
                },
            )
            ccont += 1
        plt.xlabel('Length of reads', fontsize=12)
        plt.ylabel('Density', fontsize=12)
        plt.legend(fontsize=12)
        plt.show()
        return

    def plot(self, ):
        fig, ax = plt.subplots()
        ax.bar(
            # x=[0, 1, 2, 3, 4, 5],
            # height=list(reversed([261582, 172809, 98424, 55158, 41078, 32201])),
            x=[30, 40, 50],
            height=[3, 4, 15],
            # x=[30, 40, 50],
            # height=[3579.8/60, 4579.4/60, 5553.2/60],
            width=0.2,
            color='palevioletred',
            # label=definition,
            alpha=0.75,
            # edgecolor='black',
            linewidth=0.9
        )
        ax.plot(
            # [0, 1, 2, 3, 4, 5],
            # list(reversed([261582, 172809, 98424, 55158, 41078, 32201])),
            [30, 40, 50],
            [3, 4, 15],
            # [30, 40, 50],
            # [3579.8/60, 4579.4/60, 5553.2/60],
            c='slategray',
        )
        ax.set_title('Examining the amount of reads containing potential polyA tails')
        ax.set_xlabel('Number of edit errors to be allowed in polyA tails of 200bp')
        ax.set_ylabel('Number of reads')

        # ax.set_xlabel('Number of edit errors to be allowed in polyA tails of 200bp')
        # ax.set_ylabel('Running time (minutes)')
        ax.set_xticks([30, 40, 50])
        ax.set_xticklabels([30, 40, 50], fontsize=9)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.subplots_adjust(
            # top=0.88,
            bottom=0.10,
            left=0.15,
            right=0.95,
        )
        plt.show()
        return


if __name__ == "__main__":
    DEFINE = {
        'fastq': {
            'path': to('data/'),
            # 'name': 'fastq_runid_0',
            'name': 'TSO_trimer',
        },
        'ed': [0, 1, 2],
        'std_pattern': "TTTTTTTTTTTTTTTTTTTT",

        # 'ed': [30, 40, 50],
        # 'std_pattern': 'A'*200,
    }
    p = dashboard(DEFINE)
    # print(p.pmatching())
    print(p.pmatchin1())
    # print(p.plot())