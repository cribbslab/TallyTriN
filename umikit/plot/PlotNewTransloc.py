__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None
from umikit.util.Writer import writer as gwriter
from simreadflow.util.file.read.Reader import reader as gfreader
from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from umikit.dedup.trimer.pipeline import Config
from umikit.dedup.monomer.Relation import relation as umimonorel
from simreadflow.util.random.Number import number as rannum
from Path import to


class newTransloc(Config.config):

    def __init__(self, metric, method, section='', tc_thres=None, fastq_fp=None, fn_suffix='', sv_cnt_lib_fpn=''):
        super(newTransloc, self).__init__()
        self.metric = metric
        self.tc_thres = tc_thres
        self.fn_suffix = fn_suffix
        self.gbfscc = gbfscc()
        self.gfreader = gfreader()
        self.gwriter = gwriter()
        self.rannum = rannum()
        self.umimonorel = umimonorel
        self.sv_cnt_lib_fpn = sv_cnt_lib_fpn
        self.method = method
        self.section = section
        self.fastq_fp = fastq_fp
        self.seq_num = 100
        self.num_tc_thres = 5

        if self.method == 'cnt_lib':
            with open(self.sv_cnt_lib_fpn) as fp1:
                self.cnt_dict = json.load(fp1)
            # print(self.cnt_dict)
        else:
            df = pd.DataFrame()
            # for tc_thres in np.arange(self.num_tc_thres) + 1:
            for tc_thres in [self.num_tc_thres]:
                df_real_mean, df_fake_yes_mean, df_chi_yes_mean = self.read(
                    fastq_fp=fastq_fp,
                    tc_thres=tc_thres,
                )
                df = pd.concat([
                    df,
                    df_real_mean,
                    df_fake_yes_mean,
                    df_chi_yes_mean,
                ], axis=0)
            df = df.reset_index()
            print(df)

        sns.set(font="Helvetica")
        sns.set_style("ticks")
        if section == 'plot':
            df = self.read_met_compare(fastq_fp)
            fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9, 7), sharey='all', sharex='all')
            sns.set(font="Helvetica")
            sns.set_style("ticks")

            df_noncorr = df[df['sort'] == 'noncorr']
            df_corr = df[df['sort'] == 'corr']
            tableau = [
                'royalblue',
                'purple',
                'black',
                'crimson',
            ]
            markers = [
                "o",
                "D",
                "v",
                "^",
            ]
            for i, met in enumerate(df_noncorr['method'].unique()):
                print(df_noncorr[df_noncorr['method'] == met])
                axes[0].plot(
                    np.arange(len(self.seq_errs)),
                    df_noncorr[df_noncorr['method'] == met][0],
                    label=met,
                    marker=markers[i],
                    color=tableau[i],
                    linestyle="-",
                    alpha=0.5,
                    linewidth=2
                )
                axes[1].plot(
                    np.arange(len(self.seq_errs)),
                    df_corr[df_corr['method'] == met][0],
                    label=met,
                    marker=markers[i],
                    color=tableau[i],
                    linestyle="-",
                    alpha=0.5,
                    linewidth=2
                )

            axes[0].spines['right'].set_color('none')
            axes[0].spines['top'].set_color('none')
            axes[1].spines['right'].set_color('none')
            axes[1].spines['top'].set_color('none')
            axes[0].set_ylabel('Chimeric reads (%) detected')
            axes[0].set_title('non-corrected')
            axes[1].set_ylabel('Chimeric reads (%) detected')
            axes[1].set_xlabel('Sequencing error')
            axes[1].set_title('corrected')

            axes[0].set_xticks(np.arange(len(self.seq_errs)))
            axes[1].set_xticklabels(['{:.2E}'.format(t) for t in self.seq_errs], rotation=30, ha='right', fontsize=9)
            plt.subplots_adjust(
                top=0.90,
                bottom=0.12,
                left=0.08,
                right=0.98,
                # hspace=0.40,
                # wspace=0.15,
            )
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, ncol=3, loc='upper center', bbox_transform = plt.gcf().transFigure)
            plt.show()

        if section == 'plot_line':
            # res_sns = sns.relplot(
            #     x='seq_err',
            #     y=0,
            #     hue='sort',
            #     col="tc_thres",
            #     data=df,
            #     kind="line",
            #     col_wrap=3,
            #     linewidth=2,
            #     legend="full",
            #     # estimator=np.size,
            #     # height=2.5,
            #     # ci=None,
            #     aspect=.8,
            #     # sharex=False,
            #     # palette=sns.color_palette("flare"),
            # )
            res_sns = sns.catplot(
                x='sort',
                y=0,
                # hue='sort',
                col="tc_thres",
                data=df,
                kind="violin",
                col_wrap=3,
                linewidth=0.2,
                legend="full",
                # estimator=np.size,
                # height=2.5,
                # ci=None,
                aspect=.8,
                # sharex=False,
                # palette=sns.color_palette("flare"),
            )
            # for ax in res_sns.axes.flat:
            #     ax.axvline(11, ls='--', color='tab:blue')

                # ax.axvline(11, ls='--', color='tab:blue')
                # ax.axvline(15, ls='--', color='r')
                # ax.axvline(20, ls='--', color='r')
                # ax.axvspan(15, 20, color="r", alpha=0.1)

            # leg = res_sns._legend
            # plt.setp(leg.get_title(), fontsize=14)
            # plt.setp(leg.get_texts(), fontsize=14)
            # leg.texts[0].set_text("real reads")
            # leg.texts[1].set_text("chimeric reads")
            # leg.texts[2].set_text("PCR duplicates")
            # leg.set_title("Read composition")
            # leg.set_bbox_to_anchor([0.85, 0.3])

            res_sns.set(ylim=(-0.05, 1.05))
            res_sns.set_xticklabels(rotation=30, ha='right', fontsize=9)
            res_sns.set_ylabels('Percentage of reads', fontsize=14)
            res_sns.set_xlabels('Sequencing error', fontsize=14)

            axes = res_sns.axes.flatten()
            axes[0].set_title("edit distance cutoff = 1", fontsize=16)
            axes[1].set_title("edit distance cutoff = 2", fontsize=16)
            axes[2].set_title("edit distance cutoff = 3", fontsize=16)
            axes[3].set_title("edit distance cutoff = 4", fontsize=16)
            axes[4].set_title("edit distance cutoff = 5", fontsize=16)

            plt.subplots_adjust(
                top=0.92,
                bottom=0.14,
                left=0.06,
                right=0.98,
                hspace=0.40,
                wspace=0.15,
            )
            plt.show()

        if section == 'plot_cnt':
            arr = []
            for k, v in self.cnt_dict.items():
                for k1, v1 in v.items():
                    for k2, v2 in v1.items():
                        arr.append([k, '{:.2E}'.format(float(k1)), k2, len(v2)])
            df = pd.DataFrame(arr, columns=['pn', 'seq_err', 'pos_umi', 'cnt'])
            gp = df.groupby(by=['seq_err', 'pos_umi'])
            df_gp_keys = gp.groups.keys()
            df_sub = []
            for k in df_gp_keys:
                print(k)
                df_gp = gp.get_group(k)
                if 'corr' in k[1].split('_'):
                    df_sub.append([float(k[0]), k[1], df_gp['cnt'].mean(), 'corr'])
                else:
                    df_sub.append([float(k[0]), k[1], df_gp['cnt'].mean(), 'non_corr'])
            df_sub = pd.DataFrame(df_sub, columns=['seq_err', 'pos_umi', 'cnt', 'is_corr'])
            df_sub = df_sub.sort_values(by=['seq_err'], ascending=False)
            print(df_sub)
            print(df)
            res_sns = sns.catplot(
                x='seq_err',
                y='cnt',
                hue='pos_umi',
                col="is_corr",
                data=df_sub,
                kind="bar",
                col_wrap=1,
                linewidth=0.1,
                legend="full",
                # estimator=np.size,
                # height=2.5,
                ci=None,
                aspect=1.6,
                sharey=False,
                # palette=sns.color_palette("flare"),
            )

            def change_width(ax, new_value):
                for patch in ax.patches:
                    current_width = patch.get_width()
                    diff = current_width - new_value
                    patch.set_width(new_value)
                    patch.set_x(patch.get_x() + diff * .5)

            for ax in res_sns.axes.flat:
                change_width(ax, .25)

            # leg.set_title("Read composition")
            h, l = res_sns.axes[0].get_legend_handles_labels()
            l = [
                "UMI (corrected) at 3' end",
                "merged UMI (non-corrected)",
                "merged UMI (corrected)",
                "UMI (non-corrected) at 3' end",
                "UMI (corrected) at 5' end",
                "UMI (non-corrected) at 5' end",
            ]
            res_sns._legend.remove()
            res_sns.fig.legend(h, l, ncol=3, loc='upper center')

            res_sns.set_xticklabels(["{:.2E}".format(i) for i in self.seq_errs], rotation=30, ha='right', fontsize=14)
            res_sns.set_ylabels('Number of unique UMIs', fontsize=14)
            res_sns.set_xlabels('Sequencing error', fontsize=14)

            axes = res_sns.axes.flatten()
            axes[0].set_title("corrected", fontsize=16)
            axes[1].set_title("non-corrected", fontsize=16)

            plt.subplots_adjust(
                top=0.90,
                bottom=0.14,
                left=0.06,
                right=0.98,
                hspace=0.40,
                wspace=0.15,
            )
            plt.show()

        if section == 'plot_cnt_umi':
            arr = []
            for k1, v1 in self.cnt_dict['0'].items():
                for k2, v2 in v1.items():
                    for i in v2:
                        arr.append([0, '{:.2E}'.format(float(k1)), k2, i])
            df = pd.DataFrame(arr, columns=['pn', 'sequence error', 'UMI source', 'Count'])
            print(df)
            res_sns = sns.catplot(
                x='UMI source',
                y='Count',
                # hue='UMI source',
                col="sequence error",
                data=df,
                kind="violin",
                col_wrap=4,
                linewidth=0.1,
                legend="full",
                # estimator=np.size,
                # height=2.5,
                # ci=None,
                aspect=.6,
                sharey=False,
                palette=sns.color_palette("flare"),
            )
            l = [
                "(non-corr) 3' end",
                "(non-corr) 5' end",
                "merged (non-corr)",
                "(corr) 3' end",
                "(corr) 5' end",
                "merged (corr)",
            ]
            res_sns.set_xticklabels(l, rotation=45, ha='right', fontsize=8)
            plt.subplots_adjust(
                top=0.96,
                bottom=0.12,
                left=0.06,
                right=0.98,
                hspace=0.40,
                wspace=0.15,
            )
            plt.show()

        if section == 'plot_boxplot':
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 5), sharey=False, sharex='all')
            meanpointprops = {
                'marker': 'D',
                'markeredgecolor': 'black',
                'markerfacecolor': 'black',
                'markersize': 5,
                'alpha': 0.80,
            }
            palette = [
                'pink',
                'seagreen',
            ]
            df_dc_by_cnt = self.read_boxplot(self.fastq_fp, method='dc_by_cnt', tc_thres=5)
            df_control = self.read_boxplot(self.fastq_fp, method='dc_control', tc_thres=5)
            df_cnt = df_dc_by_cnt[df_dc_by_cnt['sort'] == 'fake_yes']
            df_ctrl = df_control[df_control['sort'] == 'fake_yes']
            print(df_cnt)
            print(df_control)
            cols = [
                '10', '16',
            ]
            pos = [np.arange(len(cols))-0.2, np.arange(len(cols))+0.2]
            ctrl_coor = [-0.2, 0.2]
            for i, df_focus in enumerate([df_cnt, df_ctrl]):
                bplot = axes.boxplot(
                    df_focus[cols].values,
                    positions=pos[i],
                    showmeans=True,
                    meanprops=meanpointprops,
                    patch_artist=True,
                )
                axes.spines['right'].set_color('none')
                axes.spines['top'].set_color('none')
                axes.set_xticks([0, 1])
                axes.set_xticklabels(['{:.1E}'.format(0.005), '{:.1E}'.format(0.1)], fontsize=11)
                axes.set_xlabel('Sequencing error', fontsize=14)
                axes.set_ylabel('Chimeric reads (%) detected', fontsize=14)
                # axes.yaxis.get_major_formatter().set_powerlimits((0, 1))
                for patch in bplot['boxes']:
                    patch.set(color='black', linewidth=1)
                    # patch.set_edgecolor(color='black')
                for patch in bplot['boxes']:
                    # patch.set(facecolor='snow', alpha=0.2)
                    patch.set_facecolor('white')
                    patch.set_alpha(0.001)
                    patch.set_zorder(1)
                for whisker in bplot['whiskers']:
                    whisker.set(color='black', linewidth=2)
                for cap in bplot['caps']:
                    cap.set(color='black', linewidth=2)
                for median in bplot['medians']:
                    median.set(color='black', linewidth=2)
                for flier in bplot['fliers']:
                    flier.set(marker='o', color='y', alpha=0.5)
                for ii, i_col in enumerate(cols):
                    y = df_focus[i_col].values
                    print(len(y))
                    x = [np.random.normal(ii, 0.05, len(y)) + ctrl_coor[i], np.random.normal(ii, 0.05, len(y)) + ctrl_coor[i]]
                    axes.scatter(
                        x[ii],
                        y,
                        color=palette[i],
                        alpha=0.01,
                        s=24,
                        marker='o',
                        zorder=2,
                        # label='protein chain' if i == 0 else None,
                        # edgecolors='grey',
                        linewidths=2,
                    )
            from matplotlib.lines import Line2D
            custom_dots = [
                Line2D([0], [0], marker='o', color='w', label='umiRarity',
                          markerfacecolor='pink', markeredgecolor='pink', markersize=9, alpha=0.7),
                Line2D([0], [0], marker='o', color='w', label='control',
                       markerfacecolor='seagreen', markeredgecolor='seagreen', markersize=9, alpha=0.7),
            ]
            plt.legend(handles=custom_dots, fontsize=12)
            plt.subplots_adjust(
                top=0.92,
                bottom=0.14,
                left=0.12,
                right=0.98,
                hspace=0.40,
                wspace=0.15,
            )
            plt.show()

    def read(self, fastq_fp, tc_thres):
        df_real_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/real_yes_' + self.method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        )
        df_real_no = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/real_no_' + self.method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        )
        df_fake_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_' + self.method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        )
        df_chi_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/chi_yes_' + self.method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        )

        df_real_mean = ((df_real_yes.mean() + df_real_no.mean()) / 50).to_frame()
        df_real_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_real_mean['seq_err'] = df_real_mean.index
        df_real_mean['sort'] = 'real'
        df_real_mean['tc_thres'] = tc_thres
        print(df_real_mean)

        df_fake_yes_mean = (df_fake_yes.mean() / 250).to_frame()
        df_fake_yes_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_fake_yes_mean['seq_err'] = df_fake_yes_mean.index
        df_fake_yes_mean['sort'] = 'fake_yes'
        df_fake_yes_mean['tc_thres'] = tc_thres
        print(df_fake_yes_mean)

        df_chi_yes_mean = (df_chi_yes.mean() / 6949).to_frame()
        df_chi_yes_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_chi_yes_mean['seq_err'] = df_fake_yes_mean.index
        df_chi_yes_mean['sort'] = 'chi_yes'
        df_chi_yes_mean['tc_thres'] = tc_thres
        print(df_chi_yes_mean)

        return df_real_mean, df_fake_yes_mean, df_chi_yes_mean

    def read_boxplot(self, fastq_fp, method, tc_thres):
        df_real_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/real_yes_' + method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        ) / 250
        df_real_yes['sort'] = 'real_yes'
        df_real_no = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/real_no_' + method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        ) / 250
        df_real_no['sort'] = 'real_no'
        df_fake_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_' + method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        ) / 250
        df_fake_yes['sort'] = 'fake_yes'
        df_chi_yes = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/chi_yes_' + method + '_' + str(
                tc_thres) + '_' + self.fn_suffix + '.txt',
            header=0,
        ) / 250
        df_chi_yes['sort'] = 'chi_yes'

        df = pd.concat([
            df_real_yes,
            df_real_no,
            df_fake_yes,
            df_chi_yes,

        ], axis=0)
        df['met'] = method
        df = df.reset_index(drop=True)
        return df

    def read_met_compare(self, fastq_fp):
        df_dc_by_cnt_noncorr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_by_cnt_5_.txt',
            header=0,
        )
        df_dc_by_cnt_corr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_by_cnt_5_corr.txt',
            header=0,
        )

        df_dc_by_vote_noncorr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_by_vote_5_sm&lg.txt',
            header=0,
        )
        df_dc_by_vote_corr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_by_vote_5_sm&lg_corr.txt',
            header=0,
        )

        df_dc_control_noncorr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_control_5_.txt',
            header=0,
        )
        df_dc_control_corr = self.gfreader.generic(
            df_fpn=fastq_fp + self.metric + '/fake_yes_dc_control_5_corr.txt',
            header=0,
        )

        df_dc_by_cnt_noncorr_mean = (df_dc_by_cnt_noncorr.mean() / 250).to_frame()
        df_dc_by_cnt_noncorr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_by_cnt_noncorr_mean['seq_err'] = df_dc_by_cnt_noncorr_mean.index
        df_dc_by_cnt_noncorr_mean['method'] = 'umiRarity'
        df_dc_by_cnt_noncorr_mean['sort'] = 'noncorr'
        df_dc_by_cnt_noncorr_mean['tc_thres'] = 5
        print(df_dc_by_cnt_noncorr_mean)

        df_dc_by_cnt_corr_mean = (df_dc_by_cnt_corr.mean() / 250).to_frame()
        df_dc_by_cnt_corr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_by_cnt_corr_mean['seq_err'] = df_dc_by_cnt_corr_mean.index
        df_dc_by_cnt_corr_mean['method'] = 'umiRarity'
        df_dc_by_cnt_corr_mean['sort'] = 'corr'
        df_dc_by_cnt_corr_mean['tc_thres'] = 5
        print(df_dc_by_cnt_corr_mean)

        df_dc_by_vote_noncorr_mean = (df_dc_by_vote_noncorr.mean() / 250).to_frame()
        df_dc_by_vote_noncorr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_by_vote_noncorr_mean['seq_err'] = df_dc_by_vote_noncorr_mean.index
        df_dc_by_vote_noncorr_mean['method'] = 'umiRarityCR'
        df_dc_by_vote_noncorr_mean['sort'] = 'noncorr'
        df_dc_by_vote_noncorr_mean['tc_thres'] = 5
        print(df_dc_by_vote_noncorr_mean)

        df_dc_by_vote_corr_mean = (df_dc_by_vote_corr.mean() / 250).to_frame()
        df_dc_by_vote_corr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_by_vote_corr_mean['seq_err'] = df_dc_by_vote_corr_mean.index
        df_dc_by_vote_corr_mean['method'] = 'umiRarityCR'
        df_dc_by_vote_corr_mean['sort'] = 'corr'
        df_dc_by_vote_corr_mean['tc_thres'] = 5
        print(df_dc_by_vote_corr_mean)

        df_dc_control_noncorr_mean = (df_dc_control_noncorr.mean() / 250).to_frame()
        df_dc_control_noncorr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_control_noncorr_mean['seq_err'] = df_dc_control_noncorr_mean.index
        df_dc_control_noncorr_mean['method'] = 'control'
        df_dc_control_noncorr_mean['sort'] = 'noncorr'
        df_dc_control_noncorr_mean['tc_thres'] = 5
        print(df_dc_control_noncorr_mean)

        df_dc_control_corr_mean = (df_dc_control_corr.mean() / 250).to_frame()
        df_dc_control_corr_mean.index = ['{:.2E}'.format(i) for i in self.seq_errs]
        df_dc_control_corr_mean['seq_err'] = df_dc_control_corr_mean.index
        df_dc_control_corr_mean['method'] = 'control'
        df_dc_control_corr_mean['sort'] = 'corr'
        df_dc_control_corr_mean['tc_thres'] = 5
        print(df_dc_control_corr_mean)

        df = pd.concat([
            df_dc_by_cnt_noncorr_mean,
            df_dc_by_cnt_corr_mean,
            df_dc_by_vote_noncorr_mean,
            df_dc_by_vote_corr_mean,
            df_dc_control_noncorr_mean,
            df_dc_control_corr_mean,
        ], axis=0)
        df = df.reset_index()
        print(df)
        return df


if __name__ == "__main__":
    p = newTransloc(
        metric='seq_errs',
        method='dc_by_cnt',
        # method='dc_by_vote',
        # method='dc_control',
        # method='cnt_lib',

        # section='plot',
        # section='plot_line',
        # section='plot_cnt',
        # section='plot_cnt_umi',
        section='plot_boxplot',

        tc_thres=5,
        fn_suffix='corr',
        # fn_suffix='sm&sm',
        # fn_suffix='sm&lg',
        # fastq_fp=to('data/simu/transloc/trimer/single_read/pcr8_umi30/'),
        # sv_cnt_lib_fpn=to('data/simu/transloc/trimer/single_read/pcr8_umi30/cnt_lib.json'),

        fastq_fp=to('data/simu/transloc/trimer/single_read/1000/'),
        sv_cnt_lib_fpn=to('data/simu/transloc/trimer/single_read/1000/cnt_lib.json'),
    )