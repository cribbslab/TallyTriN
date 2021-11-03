import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umikit.util.Reader import reader as greader
from umikit.dedup.monomer.pipeline import Config
from simreadflow.gmat.FromSimulator import fromSimulator
from scipy.sparse import coo_matrix


class protSC(Config.config):

    def __init__(self, fpns):
        """
        self.gmat, _, _ = fromSimulator(simulator='SPsimSeqFixSM').run()
            self.cell_map = {k: v for k, v in enumerate(self.gmat.columns)}
            self.gene_map = {k: v for k, v in enumerate(self.gmat.index)}
            # print(self.cell_map)
            # print(self.gene_map)
            csr_ = coo_matrix(self.gmat)
            print(csr_)
            self.gbyc_arr = np.transpose([
                csr_.row.tolist(),
                csr_.col.tolist(),
                csr_.data.tolist(),
            ]).astype(np.int)
            print(pd.DataFrame(self.gbyc_arr, columns=['cell', 'gene', 'gc_cnt']).to_csv(
                to('data/simu/monomer/sc/seq_errs/ground_truth.txt'), header=True, index=False, sep='\t'))
        :param fpns:
        """
        super(protSC, self).__init__()
        self.greader = greader()
        self.df_gt = self.greader.generic(df_fpn=fpns['ground_truth'], header=0)[:]
        self.df_gt['gc_tag'] = self.df_gt.apply(lambda x: 'cell'+str(x[0] + 1) + ' ' + 'gene' + str(x[1] + 1), axis=1)
        print(self.df_gt)
        self.df = pd.DataFrame()
        df_melt = pd.DataFrame()
        for method, fpn in fpns['methods'].items():
            df_met = self.greader.generic(df_fpn=fpn, header=0)[:]
            # print(df_met)
            df_met = df_met.sub(self.df_gt['gc_cnt'], axis='rows')
            df_met = df_met.div(self.df_gt['gc_cnt'], axis='rows')
            print(df_met.div(self.df_gt['gc_cnt'], axis='rows'))


            # df_met = (df_met - 50) / 50
            df_met['gc_tag'] = self.df_gt['gc_tag']
            metric_map = {'metric' + str(k): v for k, v in enumerate(self.seq_errs)}
            df_met_melt = pd.melt(df_met, 'gc_tag', var_name='seq_err')
            df_met_melt_ = pd.DataFrame(df_met_melt, columns = ['gc_tag', 'seq_err', 'value'])
            df_met_melt_['method'] = method
            df_met_melt_['se_val'] = df_met_melt_['seq_err'].apply(lambda x: metric_map[x])
            print(df_met_melt_)
            df_met['method'] = method

            self.df = pd.concat([self.df, df_met], axis=0)
            df_melt = pd.concat([df_melt, df_met_melt_], axis=0)
        print(self.df)
        self.df_direc = self.df[self.df['method'] == 'directional']
        self.df_mcl_val = self.df[self.df['method'] == 'mcl_val']
        self.df_mcl_ed = self.df[self.df['method'] == 'mcl_ed']
        print(self.df_direc.iloc[1].values.shape)
        # grid = sns.FacetGrid(df_melt, col="gc_tag", hue="method", palette=['r', 'b', 'g'],
        #                      col_wrap=4, height=1.5)
        #
        # # Draw a horizontal line to show the starting point
        # # grid.refline(y=0, linestyle=":")
        #
        # # Draw a line plot to show the trajectory of each random walk
        # grid.map(plt.plot, "se_val", "value", marker="o")
        #
        # # Adjust the tick positions and labels
        # # grid.set(xticks=np.arange(5), yticks=[-3, 3],
        # #          xlim=(-.5, 4.5), ylim=(-3.5, 3.5))
        #
        # # Adjust the arrangement of the plots
        # grid.fig.tight_layout(w_pad=1)
        # plt.show()

    def separateline(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        row_num = 5
        col_num = 8
        fig, ax = plt.subplots(nrows=row_num, ncols=col_num, figsize=(14, 8), sharey=False, sharex='all')
        c = 0
        for i_row in range(row_num):
            for i_col in range(col_num):
                ax[i_row][i_col].plot(
                    self.seq_errs,
                    self.df_direc.iloc[c].values[:15],
                    color='tab:green',
                    alpha=0.9,
                    label='directional',
                )
                ax[i_row][i_col].plot(
                    self.seq_errs,
                    self.df_mcl_ed.iloc[c].values[:15],
                    color='tab:blue',
                    alpha=0.9,
                    label='mcl_ed',
                )
                ax[i_row][i_col].plot(
                    self.seq_errs,
                    self.df_mcl_val.iloc[c].values[:15],
                    color='crimson',
                    alpha=0.9,
                    label='mcl_val',
                )
                ax[i_row][i_col].spines['right'].set_visible(False)
                ax[i_row][i_col].spines['top'].set_visible(False)
                ax[i_row][i_col].set_title(self.df_mcl_val.iloc[c].values[15] + '\ncnt:' + str(self.df_gt.iloc[c].values[2]), fontsize=8)
                # print(self.seq_errs[::2])

                ylim_low, ylim_up = ax[i_row][i_col].get_ylim()
                # ax[i_row][i_col].set_ylim(ylim_low, ylim_up + 0.01)
                # print(ylim_low)
                ylims = [round(i) for i in np.linspace(0, ylim_up, 3)]
                ax[i_row][i_col].set_yticks(ylims)
                ax[i_row][i_col].set_yticklabels(ylims, fontsize=8)

                xlims = np.linspace(self.seq_errs[0], self.seq_errs[-1], 4)
                ax[i_row][i_col].set_xticks(xlims)
                ax[i_row][i_col].set_xticklabels(['{:.1e}'.format(i) for i in xlims], rotation=20, fontsize=8)
                # ax[i_row][i_col].legend(loc='upper left', fontsize=5)
                c += 1

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Sequencing error rate', fontsize=12)
        plt.ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)

        # fig.supxlabel('Sequencing error rate', fontsize=12)
        # fig.supylabel(r'$\frac{N_e-N_t}{N_t}$')
        # plt.tight_layout()
        fig.subplots_adjust(
            # top=0.98,
            # bottom=0.15,
            # left=0.15,
            # right=0.98,
            hspace=0.70,
            wspace=0.35
        )

        lines = []
        labels = []
        axLine, axLabel = ax[0][0].get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
        plt.legend(lines, labels, fontsize=9, bbox_to_anchor=(0.5, 1.15), loc=9, ncol=7)
        plt.show()
        return


if __name__ == "__main__":
    from Path import to
    DEFINE = {
        'fpns': {
            'methods': {
                'directional': to('data/simu/monomer/sc/seq_errs/directional.txt'),
                'mcl_ed': to('data/simu/monomer/sc/seq_errs/mcl_ed.txt'),
                'mcl_val': to('data/simu/monomer/sc/seq_errs/mcl_val.txt'),
            },
            'ground_truth': to('data/simu/monomer/sc/seq_errs/ground_truth.txt'),
        },


    }
    p = protSC(DEFINE['fpns'])

    print(p.separateline())