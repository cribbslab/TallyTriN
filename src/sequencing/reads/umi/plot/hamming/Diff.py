__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from src.sequencing.reads.umi.plot.hamming.Style import style
from src.util.file.read.Reader import reader as gfreader
from src.util.file.write.Writer import writer as fwriter
from Path import to
from matplotlib import rcParams


class diff(object):

    def __init__(self, ):
        self.gfreader = gfreader()
        self.fwriter = fwriter()
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Tahoma']
        self.df_monomer = self.gfreader.generic(
            df_fpn=to('data/simu/umi/seq_errs/monomer/trimmed/seq.txt')
        )
        self.df_dimer = self.gfreader.generic(
            df_fpn=to('data/simu/umi/seq_errs/dimer/trimmed/seq.txt')
        )
        self.df_trimer = self.gfreader.generic(
            df_fpn=to('data/simu/umi/seq_errs/trimer/trimmed/seq.txt')
        )
        print(self.df_monomer.loc[:, 14:].values.T)
        print(self.df_dimer.loc[:, 14:].values.T)
        print(self.df_trimer.loc[:, 14:].values.T)
        self.met_val_dict = {
            'monomer': self.df_monomer.loc[:, 11:].values.T,
            'dimer': self.df_dimer.loc[:, 11:].values.T,
            'trimer': self.df_trimer.loc[:, 11:].values.T,
        }
        self.met_pos_dict = {}
        for i, (k, v) in enumerate(self.met_val_dict.items()):
            self.met_pos_dict[k] = np.arange((19-11)) * 3 - i * 0.6

    def boxplot(self, ):
        meanpointprops = {
            'marker': 'o',
            'markeredgecolor': 'grey',
            'markerfacecolor': 'grey',
            'markersize': 6
        }

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11, 5))
        bplot_handles = []        # bplot_handles store methods with size equal to the number of the methods
        for val_key, pos_key in zip(self.met_val_dict, self.met_pos_dict):
            print(val_key)
            val_arr_2d = self.met_val_dict[val_key]
            print(val_arr_2d.shape)
            pos_arr_1d = self.met_pos_dict[pos_key]
            # print(pos_arr_1d)
            bplot_handles.append(ax.boxplot(
            x=list(val_arr_2d),
            positions=list(pos_arr_1d),
            showmeans=True,
            meanprops=meanpointprops,
            patch_artist=True,
            showfliers=False, # hide outliers, the same as sym=''
            ))
        x_ticks = np.arange(0, (19-11) * 3, 3)
        tt = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
              0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3]
        x_tick_labels = ["{:.1e}".format(i) for i in tt][11:]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, rotation=30, ha='right', fontsize=9)

        from matplotlib import rcParams
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Tahoma']
        ax.set_xlabel('Sequencing error rate', fontsize=18)
        ax.set_ylabel('Number of incorrect UMIs', fontsize=18)
        # ax.set_ylabel('Number of nucleotide errors', fontsize=18)
        ax.yaxis.grid(True)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        fig.subplots_adjust(
            bottom=0.16,
            # top=0.92, left=0.10, right=0.95, hspace=0.25, wspace=0.17
        )
        palette = [
            'magenta',
            'skyblue',
            'gray',
            'sienna',
            'crimson',
        ]
        style().plain(boxplot_handles=bplot_handles, palette=palette)

        # add legend to the boxplot
        ax.legend(
            handles=[i["boxes"][0] for i in bplot_handles],
            labels=[*self.met_val_dict.keys()],
            bbox_to_anchor=(0.5, 1.15),
            loc=9,
            ncol=7,
        )
        # fig.tight_layout()
        plt.show()
        return 0

    def linebase(self, ):
        stf1 = self.minmax(self.df_monomer.sum(axis=0))
        stf2 = self.minmax(self.df_dimer.sum(axis=0))
        stf3 = self.minmax(self.df_trimer.sum(axis=0))
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
        x_ticks = np.arange(19)
        ax.plot(x_ticks, stf1, label='monomer', c='orange', lw=3, alpha=0.7)
        ax.plot(x_ticks, stf2, label='dimer', c='royalblue', lw=3, alpha=0.7)
        ax.plot(x_ticks, stf3, label='trimer', c='firebrick', lw=3, alpha=0.7)
        tt = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
         0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3]
        x_tick_labels = ["{:.1e}".format(i) for i in tt]
        ax.set_xlabel('Sequencing error rate', fontsize=13)
        ax.set_ylabel('Number of nucleotide errors', fontsize=13)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, rotation=30, ha='right', fontsize=9)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend()

        fig.subplots_adjust(
            top=0.92,
            bottom=0.18,
            left=0.11,
            right=0.95,
            hspace=0.20,
            wspace=0.15
        )
        plt.show()
        return

    def svseq(self, ):
        cols = self.df_monomer.columns
        for i in cols:
            print(i)
            self.df_monomer[i] = self.df_monomer[i].apply(lambda x: 1 if x != 0 else 0)
            self.df_dimer[i] = self.df_dimer[i].apply(lambda x: 1 if x != 0 else 0)
            self.df_trimer[i] = self.df_trimer[i].apply(lambda x: 1 if x != 0 else 0)
        print(self.df_monomer)
        print(self.df_dimer)
        print(self.df_trimer)
        self.fwriter.generic(df=self.df_monomer, sv_fpn=to('data/simu/umi/seq_errs/monomer/trimmed/seq.txt'), df_sep='\t')
        self.fwriter.generic(df=self.df_dimer, sv_fpn=to('data/simu/umi/seq_errs/dimer/trimmed/seq.txt'), df_sep='\t')
        self.fwriter.generic(df=self.df_trimer, sv_fpn=to('data/simu/umi/seq_errs/trimer/trimmed/seq.txt'), df_sep='\t')
        return

    def lineseq(self, ):
        stf1 = self.minmax(self.df_monomer.sum(axis=0))
        stf2 = self.minmax(self.df_dimer.sum(axis=0))
        stf3 = self.minmax(self.df_trimer.sum(axis=0))
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
        x_ticks = np.arange(19)
        ax.plot(x_ticks, stf1, label='monomer', c='orange', lw=3, alpha=0.7)
        ax.plot(x_ticks, stf2, label='dimer', c='royalblue', lw=3, alpha=0.7)
        ax.plot(x_ticks, stf3, label='trimer', c='firebrick', lw=3, alpha=0.7)
        tt = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
              0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3]
        x_tick_labels = ["{:.1e}".format(i) for i in tt]
        ax.set_xlabel('Sequencing error rate', fontsize=13)
        ax.set_ylabel('Number of incorrect UMIs', fontsize=13)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, rotation=30, ha='right', fontsize=9)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend()

        fig.subplots_adjust(
            top=0.92,
            bottom=0.18,
            left=0.11,
            right=0.95,
            hspace=0.20,
            wspace=0.15
        )
        plt.show()
        return

    def minmax(self, arr):
        return arr

    def minmax1(self, arr):
        max_ = max(arr)
        return [i/max_ for i in arr]

    def stackedbar(self, ):
        tt = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
              0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3]
        x_tick_labels = ["{:.1e}".format(i) for i in tt]
        index = pd.Index(x_tick_labels, name='test')

        stf1 = self.minmax(self.df_monomer.sum(axis=0))
        stf2 = self.minmax(self.df_dimer.sum(axis=0))
        stf3 = self.minmax(self.df_trimer.sum(axis=0))

        data = {'monomer': stf1.values.tolist(),
                'dimer': stf2.values.tolist(),
                'trimer': stf3.values.tolist(),
                }

        df = pd.DataFrame(data, index=index)
        ax = df.plot(
            kind='bar', stacked=True, figsize=(10, 6),
            color={
            "monomer": "magenta", "dimer": "skyblue", "trimer": "gray",
            },
            # color='None',
            # colormap='Paired',
            alpha=0.4,
        )
        ax.set_xlabel('Sequencing error rate', fontsize=15)
        ax.set_ylabel('Number of incorrect UMIs', fontsize=15)
        # ax.set_ylabel('Number of nucleotide errors', fontsize=15)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.subplots_adjust(
            top=0.92,
            bottom=0.18,
            left=0.11,
            right=0.85,
            hspace=0.20,
            wspace=0.15
        )
        plt.legend(title='labels', bbox_to_anchor=(1.0, 1), loc='upper left')
        # plt.savefig('stacked.png')  # if needed
        plt.show()
        return


if __name__ == "__main__":

    p = diff()
    # print(p.boxplot())
    # print(p.linebase())
    # print(p.lineseq())
    print(p.stackedbar())