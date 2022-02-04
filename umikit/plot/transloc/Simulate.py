import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umikit.util.Reader import reader as greader
from umikit.dedup.trimer.pipeline import Config


class protTrimerSplit(Config.config):

    def __init__(self, fpns):
        super(protTrimerSplit, self).__init__()
        self.greader = greader()
        self.df = pd.DataFrame()
        self.df_T = pd.DataFrame()
        for method, fpn in fpns.items():
            df_met = self.greader.generic(df_fpn=fpn, header=0)[:]
            df_met_T = df_met.T
            df_met_T = (df_met_T - 50) / 50
            df_met_T.columns = ['{:.1e}'.format(x) for x in self.seq_errs[:]]
            # df_met_T.columns = [int(x*100000) for x in self.seq_errs[:]]
            df_met_T['method'] = method
            # print(df_met_T)
            # df_met = np.exp((df_met - 50) / 50)
            # df_met = (df_met - 50) / 50
            df_met['mean'] = df_met.mean(axis=1)
            df_met['max'] = df_met.max(axis=1)
            df_met['min'] = df_met.min(axis=1)
            df_met['std'] = df_met.std(axis=1)
            df_met['mean-min'] = df_met['std']
            df_met['max-mean'] = df_met['std']
            df_met['method'] = method
            df_met['metric'] = ['{:.1e}'.format(x) for x in self.seq_errs[:]]
            # df_met['metric'] = ['{:.1f}'.format(x) for x in self.ampl_rates[:]]
            # df_met['metric'] = ['{:.0f}'.format(x) for x in self.umi_unit_lens[:]]
            # df_met['metric'] = ['{:.1e}'.format(x) for x in self.pcr_errs[:]]
            self.df = pd.concat([self.df, df_met], axis=0)
            self.df_T = pd.concat([self.df_T, df_met_T], axis=0)
        print(self.df)
        self.df_direc = self.df[self.df['method'] == 'directional']
        self.df_mcl_val = self.df[self.df['method'] == 'mcl_val']
        self.df_mcl_ed = self.df[self.df['method'] == 'mcl_ed']
        # print(self.df_T)
        # self.df = self.df.reset_index(drop=True)

        # sns.jointplot(data=self.df, x="mean-min", y="max", hue="met",)
        plt.show()
        self.df_melt = pd.melt(self.df_T, 'method', var_name="Sequencing error")
        # self.df_melt = pd.DataFrame(self.df_melt)
        print(self.df_melt)
        self.df_gp = self.df.groupby(by=['method'])
        self.gp_keys = self.df_gp.groups.keys()
        # print(self.gp_keys)

    def jointplot(self, ):
        # fig, ax = plt.subplots()
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        ppp = sns.jointplot(
            # x=self.df_mcl_ed['mean'].values,
            x=self.df_mcl_val['mean'].values,
            y=self.df_direc['mean'].values,
            kind="reg",
            color="crimson",
            label='asd',
        )
        ppp.ax_joint.plot([0, 50], [0, 50], 'grey', linewidth=2, alpha=1)
        # ppp.set_axis_labels('mcl_ed '+ r'($\frac{N_e-N_t}{N_t}$)' , 'directional ' + r'($\frac{N_e-N_t}{N_t}$)', fontsize=14)
        ppp.set_axis_labels('mcl_val '+ r'($\frac{N_e-N_t}{N_t}$)' , 'directional ' + r'($\frac{N_e-N_t}{N_t}$)', fontsize=14)
        ppp.ax_joint.text(40, 45, "confidence interval", horizontalalignment='right', size='medium', color='crimson',)
        ppp.ax_joint.text(20, 35, "regression", horizontalalignment='left', size='medium', color='crimson',)
        ppp.ax_joint.text(34, 32, "baseline", horizontalalignment='left', size='medium', color='black',)
        sns.despine(right=True, top=True)
        plt.tight_layout()
        plt.show()

    def errorbar(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        # plt.rc('text', usetex=True)
        # plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
        fig, ax = plt.subplots()
        # ax.plot(
        #     file_n301.index,
        #     n301_means,
        #     label='Full TrainData dataset',
        #     alpha=0.9,
        #     linewidth=3.0,
        #     # s=,
        #     c='royalblue'
        # )
        cc = [
            'tab:green',
            'tab:blue',
            'crimson',
        ]

        for i, met in enumerate(self.gp_keys):
            df_met = self.df_gp.get_group(met)

            ax.errorbar(
                x=df_met.index,
                y=df_met['mean'],
                yerr=[df_met['mean-min'], df_met['max-mean']],
                fmt='o',
                alpha=0.6,
                ecolor=cc[i],
                color=cc[i],
                linewidth=1,
                elinewidth=1,
                capsize=2,
                markersize=6,
                label=met,
            )
            # print(df_met['metric'])
            ax.set_xticks(np.arange(df_met['metric'].shape[0]))
            ax.set_xticklabels(df_met['metric'], fontsize=8, rotation=45)
        # ax.set_xlabel('Sequencing error rate', fontsize=12)
        # ax.set_xlabel('Amplification rate', fontsize=12)
        # ax.set_xlabel('UMI length', fontsize=12)
        ax.set_xlabel('Polymerase error rate', fontsize=12)
        ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=16)
        # ax.set_title(DEFINE['title'], fontsize=12)
        sns.despine(right=True, top=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.16,
            left=0.11,
            # left=0.12,
            right=0.95,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.legend(fontsize=11, loc='upper left')
        plt.show()
        return

    def errorband(self, ):
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        cc = [
            # 'tab:green',
            # 'tab:blue',
            # 'crimson',
            # 'tab:orange',
            # 'tab:brown',

            'midnightblue',
            'cornflowerblue',
            'maroon',
            'lightcoral',
        ]
        for i, met in enumerate(self.gp_keys):
            df_met = self.df_gp.get_group(met)
            ax.errorbar(
                x=df_met.index,
                y=df_met['mean'],
                yerr=[df_met['mean-min'], df_met['max-mean']],
                fmt='o',
                alpha=0.8,
                ecolor=cc[i],
                color=cc[i],
                linestyle='-',
                linewidth=2,
                elinewidth=0.5,
                capsize=2,
                markersize=3,
                label=met,
            )
            ax.plot(df_met.index, df_met['mean'] - df_met['mean-min'], color=cc[i], linewidth=0.1, alpha=0.1)
            ax.plot(df_met.index, df_met['max-mean'] + df_met['mean'], color=cc[i], linewidth=0.1, alpha=0.1)
            ax.fill_between(
                df_met.index,
                df_met['mean'] - df_met['mean-min'],
                df_met['max-mean'] + df_met['mean'],
                alpha=0.1,
                color=cc[i],
            )
            ax.set_xticks(df_met['metric'].index)
            ax.set_xticklabels(df_met['metric'], fontsize=8, rotation=45)

        ax.set_xlabel('Sequencing error rate', fontsize=12)
        ax.set_ylim(-3.0, )
        # ax.set_xlabel('Amplification rate', fontsize=12)
        ax.set_ylabel('Number of molecules failing to be \n estimated as artifactual/real', fontsize=12)
        # ax.set_ylabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=14)
        # ax.set_title(DEFINE['title'], fontsize=12)
        sns.despine(right=True, top=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.16,
            left=0.13,
            right=0.95,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.legend(fontsize=11)
        plt.show()
        return


if __name__ == "__main__":
    from Path import to
    # metric_char ='pcr_nums'
    # metric_char ='pcr_errs'
    metric_char ='seq_errs'
    # metric_char  = 'ampl_rates'
    # metric_char ='umi_lens'

    DEFINE = {
        'fpns': {
            # 'artifactual (actual)': to('data/simu/transloc/trimer/single_read/') + metric_char + '/act_fake_healing_sgl.txt',
            # 'artifactual (estimated)': to('data/simu/transloc/trimer/single_read/') + metric_char + '/est_fake_healing_sgl.txt',
            #
            # 'real (actual)': to('data/simu/transloc/trimer/single_read/') + metric_char + '/act_real_healing_sgl.txt',
            # 'real (estimated)': to('data/simu/transloc/trimer/single_read/') + metric_char + '/est_real_healing_sgl.txt',

            'artifactual (actual)': to('data/simu/transloc/trimer/bulk/') + metric_char + '/act_fake_healing_bulk.txt',
            'artifactual (estimated)': to('data/simu/transloc/trimer/bulk/') + metric_char + '/est_fake_healing_bulk.txt',

            'real (actual)': to('data/simu/transloc/trimer/bulk/') + metric_char + '/act_real_healing_bulk.txt',
            'real (estimated)': to('data/simu/transloc/trimer/bulk/') + metric_char + '/est_real_healing_bulk.txt',
        },
    }
    p = protTrimerSplit(
        fpns=DEFINE['fpns']
    )
    # print(p.errorbar())
    print(p.errorband())