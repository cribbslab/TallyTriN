import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umikit.util.Reader import reader as greader
from umikit.dedup.monomer.pipeline import Config


class protocol(Config.config):

    def __init__(self, fpns):
        super(protocol, self).__init__()
        self.greader = greader()
        self.df = pd.DataFrame()
        self.df_T = pd.DataFrame()
        for method, fpn in fpns.items():
            df_met = self.greader.generic(df_fpn=fpn, header=0)[:]
            df_met_T = df_met.T
            df_met_T = (df_met_T - 50) / 50
            # df_met_T.columns = ['{:.1e}'.format(x) for x in self.seq_fix_errs8:9]]
            df_met_T.columns = [int(x*100000) for x in self.seq_fix_errs[:]]
            df_met_T['met'] = method
            # print(df_met_T)
            df_met = (df_met - 50) / 50
            df_met['mean'] = df_met.mean(axis=1)
            df_met['max'] = df_met.max(axis=1)
            df_met['min'] = df_met.min(axis=1)
            df_met['mean-min'] = df_met['mean'] - df_met['min']
            df_met['max-mean'] = df_met['max'] - df_met['mean']
            df_met['min'] = df_met.min(axis=1)
            df_met['std'] = df_met.std(axis=1)
            df_met['met'] = method
            df_met['metric'] = ['{:.1e}'.format(x) for x in self.seq_fix_errs[:]]
            self.df = pd.concat([self.df, df_met], axis=0)
            self.df_T = pd.concat([self.df_T, df_met_T], axis=0)
        print(self.df)
        self.df_direc = self.df[self.df['met'] == 'directional']
        self.df_mcl_val = self.df[self.df['met'] == 'mcl_val']
        self.df_mcl_ed = self.df[self.df['met'] == 'mcl_ed']
        # print(self.df_T)
        self.df = self.df.reset_index(drop=True)

        # sns.jointplot(data=self.df, x="mean-min", y="max", hue="met",)
        plt.show()
        self.df_melt = pd.melt(self.df_T, 'met', var_name="measurement")
        # self.df_melt = pd.DataFrame(self.df_melt)
        print(self.df_melt)
        self.df_gp = self.df.groupby(by=['met'])
        self.gp_keys = self.df_gp.groups.keys()
        # print(self.gp_keys)

    def jointplot(self, ):
        # fig, ax = plt.subplots()
        sns.set(font="Helvetica")
        sns.set_style("ticks")
        ppp = sns.jointplot(
            x=self.df_mcl_ed['mean'].values,
            y=self.df_direc['mean'].values,
            kind="reg",
            color="crimson",
            label='asd',
        )
        ppp.ax_joint.plot([0, 50], [0, 50], 'grey', linewidth=2, alpha=1)
        ppp.set_axis_labels('mcl_ed', 'directional', fontsize=12)
        sns.despine(right=True, top=True)
        plt.tight_layout()
        plt.show()

    def jointgrid(self, ):
        ddd = pd.DataFrame()
        ddd['met'] = self.df_melt['met']
        ddd['measurement'] = self.df_melt['measurement']
        ddd['value'] = self.df_melt['value']
        print(self.df_melt['measurement'].unique())
        # g = sns.JointGrid(data=ddd[ddd['met'] == 'mcl_ed'], x="measurement", y="value", marginal_ticks=True)
        # g = sns.JointGrid(data=ddd[ddd['met'] == 'mcl_val'], x="measurement", y="value", marginal_ticks=True)
        g = sns.JointGrid(data=ddd[ddd['met'] == 'directional'], x="measurement", y="value", marginal_ticks=True)

        # Set a log scaling on the y axis
        # g.ax_joint.set(yscale="log")
        # Create an inset legend for the histogram colorbar
        cax = g.figure.add_axes([.15, .55, .02, .2])

        # Add the joint and marginal histogram plots
        g.plot_joint(
            sns.kdeplot, discrete=(True, False),
            cmap="light:#03012d", pmax=.8, cbar=True, cbar_ax=cax
        )
        g.plot_marginals(
            sns.histplot,
            element="step",
            color="#03012d",
        )
        # ax.legend(ncol=2, loc="upper right", frameon=True)
        # # ax.set(ylabel="", xlabel="Automobile collisions per billion miles")
        # ax.set_ylabel('Sequencing error', fontsize=12)
        # ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=12)
        # sns.despine(left=True, bottom=True)
        # fig.subplots_adjust(
        #     top=0.98,
        #     bottom=0.12,
        #     left=0.20,
        #     right=0.98,
        #     # hspace=0.40,
        #     # wspace=0.15
        # )
        plt.show()
        return

    def strip(self, ):
        # f, ax = plt.subplots()
        # sns.despine(bottom=True, left=True)
        #
        # # Show each observation with a scatterplot
        # sns.stripplot(x="value", y="measurement", hue="met",
        #               data=self.df_melt, dodge=True, alpha=.25, zorder=1)
        # sns.pointplot(x="value", y="measurement", hue="met",
        #               data=self.df_melt, dodge=.8 - .8 / 3,
        #               join=False, palette="dark",
        #               markers="d", scale=.75, ci=None)
        #
        # # Improve the legend
        # handles, labels = ax.get_legend_handles_labels()
        # ax.legend(handles[3:], labels[3:], title="met",
        #           handletextpad=0, columnspacing=1,
        #           loc="upper right", ncol=3, frameon=True)
        plt.show()

    def stackedbar(self, ):
        fig, ax = plt.subplots(figsize=(4, 5))
        # sns.set(font="Verdana")
        sns.set(font="Helvetica")

        self.df_mcl_val["dmean"] = self.df_direc["mean"] - self.df_mcl_val["mean"]
        self.df_mcl_ed["dmean"] = self.df_direc["mean"] - self.df_mcl_ed["mean"]

        # self.df_mcl_val["dmean"] = np.exp(self.df_direc["mean"] - self.df_mcl_val["mean"])
        # self.df_mcl_ed["dmean"] = np.exp(self.df_direc["mean"] - self.df_mcl_ed["mean"])

        sns.set_color_codes("pastel")
        sns.barplot(
            x="dmean",
            y="metric",
            data=self.df_mcl_val,
            label="dFC_ed",
            color="b",
        )
        sns.set_color_codes("muted")
        sns.barplot(
            x="dmean",
            y="metric",
            data=self.df_mcl_ed,
            label="dFC_val",
            color="b",
        )

        # plt.fill_between(self.df_mcl_val["dmean"], self.df_mcl_val["metric"])
        ax.legend(ncol=2, loc="upper right", frameon=True)
        # ax.set(ylabel="", xlabel="Automobile collisions per billion miles")
        ax.set_ylabel('Sequencing error', fontsize=12)
        ax.set_xlabel(r'$\frac{N_e-N_t}{N_t}$', fontsize=12)
        sns.despine(left=True, bottom=True)
        fig.subplots_adjust(
            top=0.98,
            bottom=0.12,
            left=0.20,
            right=0.98,
            # hspace=0.40,
            # wspace=0.15
        )
        plt.show()

    def errorbar(self, ):
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
                label=met,
            )
        # ax.set_xticks(np.arange(0, 100, 10))
        # ax.set_xticklabels(np.arange(1, 101, 10), fontsize=9)
        # ax.set_xlabel(DEFINE['xlabel'], fontsize=10)
        # ax.set_ylabel(DEFINE['ylabel'], fontsize=10)
        # ax.set_title(DEFINE['title'], fontsize=12)
        plt.legend(fontsize=11)
        plt.show()
        return

    def errorband(self, ):
        fig, ax = plt.subplots()
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
                # fmt='o',
                alpha=0.7,
                ecolor=cc[i],
                color=cc[i],
                linewidth=2,
                elinewidth=0.5,
                capsize=2,
                label=met,
            )
            ax.plot(df_met.index, df_met['mean'] - df_met['mean-min'], color=cc[i], alpha=0.05)
            ax.plot(df_met.index, df_met['max-mean'] + df_met['mean'], color=cc[i], alpha=0.05)
            ax.fill_between(
                df_met.index,
                df_met['mean'] - df_met['mean-min'],
                df_met['max-mean'] + df_met['mean'],
                alpha=0.1,
                color=cc[i],
            )
        # ax.set_xticks(np.arange(0, 100, 10))
        # ax.set_xticklabels(np.arange(1, 101, 10), fontsize=9)
        # ax.set_xlabel(DEFINE['xlabel'], fontsize=10)
        # ax.set_ylabel(DEFINE['ylabel'], fontsize=10)
        # ax.set_title(DEFINE['title'], fontsize=12)
        plt.legend(fontsize=11)
        plt.show()
        return


if __name__ == "__main__":
    from Path import to
    DEFINE = {
        'fpns': {
            'directional': to('data/simu/monomer/general/ar1/seq_errs/directional.txt'),
            'mcl_val': to('data/simu/monomer/general/ar1/seq_errs/mcl_val.txt'),
            'mcl_ed': to('data/simu/monomer/general/ar1/seq_errs/mcl_ed.txt'),
        },
    }
    p = protocol(
        fpns=DEFINE['fpns']
    )
    print(p.jointplot())
    print(p.jointgrid())
    # print(p.stackedbar())
    # print(p.errorbar())
    # print(p.errorband())