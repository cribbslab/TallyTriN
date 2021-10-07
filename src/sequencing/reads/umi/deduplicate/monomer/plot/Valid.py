import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx


class valid(object):

    def __init__(self, ):
        pass

    def n1(self, df_disapv, df_apv):
        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax[0].plot(
            df_disapv.index,
            df_disapv['0'],
            label='different',
            color='cornflowerblue',
            lw=3,
        )
        ax[0].plot(
            df_disapv.index,
            df_disapv['1'],
            label='same',
            color='salmon',
            lw=3,
        )
        ax[1].plot(
            df_apv.index,
            df_apv['0'],
            label='different',
            color='cornflowerblue',
            lw=3,
        )
        ax[1].plot(
            df_apv.index,
            df_apv['1'],
            label='same',
            color='salmon',
            lw=3,
        )

        # ax[0].set_xlabel('Time (ps)', fontsize=14)

        ax[0].set_ylabel('number of pairwise \nunique UMIs', fontsize=11)
        ax[0].set_title('Not merged', fontsize=12)
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['top'].set_visible(False)

        # ax[1].set_xlabel('Time (ps)', fontsize=14)
        ax[1].set_xticks(df_apv.index)
        ax[1].set_xticklabels(df_apv['metric'].apply(lambda x: 'PCR #' + x), fontsize=7, rotation=30)
        ax[1].set_ylabel('number of pairwise \nunique UMIs', fontsize=11)
        ax[1].set_title('Merged', fontsize=12)
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        # sns.lineplot(data=data, palette="tab10", linewidth=2.5)
        handles1, labels1 = ax[0].get_legend_handles_labels()
        ax[0].legend(
            handles1,
            labels1,
            fontsize=10,
        )

        handles2, labels2 = ax[1].get_legend_handles_labels()
        ax[1].legend(
            handles2,
            labels2,
            fontsize=10,
        )

        fig.subplots_adjust(
            # top=0.92,
            # bottom=0.13,
            # left=0.13,
            # right=0.95,
            hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def n2(self, df):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
        df_gp = df.groupby(by=['method'])
        df_gp_keys = df_gp.groups.keys()
        method_map = {
            'ccs': 'cluster',
            'adj': 'adjacency',
            'direc': r'$directional$',
            'mcl_val': 'MCL-val',
            'mcl_ed': 'MCL-ed',
        }
        palette = {
            'ccs': 'red',
            'adj': 'black',
            'direc': 'steelblue',
            'mcl_val': 'chocolate',
            'mcl_ed': 'firebrick',
            # 'black',
            # 'chocolate',
            # 'saddlebrown',
            # 'darkgoldenrod',
            # 'firebrick',
        }
        for method in df_gp_keys:
            df_met = df_gp.get_group(method)
            print(df_met)
            # if method != 'adj':
            ax.plot(
                df_met['metric'],
                df_met['dedup_cnt'].apply(
                    # lambda x: x / 50
                    lambda x: (x - 50) / 50
                    # lambda x: np.exp((x - 50) / 50)
                ),
                label=method_map[method],
                color=palette[method],
                lw=2.5,
                alpha=0.7
            )
        # ax[0].set_xlabel('Time (ps)', fontsize=14)
        c = df.loc[df['method'] == 'mcl_ed']['metric']
        ax.set_xticks(c)
        ax.set_xticklabels(c, fontsize=8)
        # ax.set_xticklabels(c.astype(np.float).apply(lambda x: '{:.2e}'.format(x)), fontsize=8, rotation=30)
        # ax.set_xticklabels(c.astype(np.float).round(1), fontsize=8)

        # ax.set_xlabel('PCR cycle', fontsize=11)
        # ax.set_xlabel('Polymerase error', fontsize=11)
        # ax.set_xlabel('Sequencing error', fontsize=11)
        ax.set_xlabel('UMI length', fontsize=11)
        # ax.set_xlabel('Amplification rate', fontsize=11)
        ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # sns.lineplot(data=data, palette="tab10", linewidth=2.5)
        handles1, labels1 = ax.get_legend_handles_labels()
        ax.legend(
            handles1,
            labels1,
            fontsize=10,
        )
        fig.subplots_adjust(
            # top=0.92,
            # bottom=0.13,
            # left=0.13,
            # right=0.95,
            hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def n2dist(self, df):
        sns.displot(data=df, x='dedup_cnt', hue='method', kind="kde", rug=True)
        plt.show()