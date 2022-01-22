__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import pandas as pd
from umikit.fastq.Convert import convert as fas2bam
from umikit.trim.Template import template as umitrim
from umikit.util.Writer import writer as gwriter
from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from umikit.dedup.dimer.pipeline import Config
from umikit.dedup.monomer.Relation import relation as umimonorel
from umikit.deduplicate.monomer.DedupPos import dedupPos
from umikit.plot.Valid import valid as plotv
from Path import to


class umi(Config.config):

    def __init__(self, metric, method, comp_cat, fastq_fp=None, is_trim=False, is_tobam=False, is_dedup=False):
        super(umi, self).__init__()
        self.metric = metric
        self.comp_cat = comp_cat
        self.gbfscc = gbfscc()
        self.umimonorel = umimonorel
        self.method = method
        self.gwriter = gwriter()
        self.plotv = plotv()
        df_dedup = pd.DataFrame()
        for i_pn in range(self.permutation_num):
            dedup_arr = []
            for id, i_metric in enumerate(self.metric_vals[self.metric]):
                if self.metric == 'pcr_nums':
                    print('=>at PCR {}'.format(i_metric))
                    fn_surf = str(i_metric)
                    self.mcl_inflat = 2.3
                    self.mcl_exp = 2
                    self.mcl_fold_thres = 1
                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                elif self.metric == 'pcr_errs':
                    self.mcl_inflat = i_metric
                    print('=>No.{} PCR error: {}'.format(id, i_metric))
                    fn_surf = str(id)
                    # # /*** mcl_ed params ***/
                    # self.mcl_inflat = 1.1 if i_metric > 0.005 else 1.7
                    # self.mcl_exp = 2
                    # self.mcl_fold_thres = 1

                    # # /*** mcl_val params ***/
                    self.mcl_inflat = 1.1 if i_metric > 0.005 else 1.8
                    self.mcl_exp = 2
                    self.mcl_fold_thres = 2
                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                elif self.metric == 'seq_errs':
                    print('=>No.{} sequencing error: {}'.format(id, i_metric))
                    # self.mcl_inflat = 1.1 if i_metric > 0.005 else 2.7
                    # self.mcl_exp = 3
                    self.mcl_fold_thres=1.6
                    self.mcl_inflat = 1.1 if i_metric > 0.005 else 2.7
                    self.mcl_exp = 2
                    fn_surf = str(id)
                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                elif self.metric == 'ampl_rates':
                    print('=>No.{} amplification rate: {}'.format(id, i_metric))
                    fn_surf = str(id)
                    # self.mcl_inflat = 1.3 if i_metric > 0.5 else 2
                    # # /*** mcl_ed params ***/
                    # if i_metric < 8:
                    #     self.mcl_inflat = 4
                    # if i_metric >= 8 and i_metric <= 11:
                    #     self.mcl_inflat = 2.3
                    # if i_metric > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 3
                    # self.mcl_fold_thres = 1

                    # /*** mcl_val params ***/
                    if i_metric < 8:
                        self.mcl_inflat = 2
                    if i_metric >= 0.9:
                        self.mcl_inflat = 1.8
                    self.mcl_exp = 4
                    self.mcl_fold_thres = 11

                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                elif self.metric == 'umi_lens':
                    print('=>No.{} UMI length: {}'.format(id, i_metric))
                    fn_surf = str(i_metric)
                    # self.mcl_inflat = 1.1 if i_metric > 11 else 2.3
                    # # /*** mcl_ed params ***/
                    # if i_metric < 8:
                    #     self.mcl_inflat = 4
                    # if i_metric >= 8 and i_metric <= 11:
                    #     self.mcl_inflat = 2.3
                    # if i_metric > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 3
                    # self.mcl_fold_thres = 1

                    # # # /*** mcl_val params ***/
                    # if i_metric < 8:
                    #     self.mcl_inflat = 6
                    # if i_metric >= 8 and i_metric <= 11:
                    #     self.mcl_inflat = 2.3
                    # if i_metric > 11:
                    #     self.mcl_inflat = 1.1
                    # self.mcl_exp = 4
                    # self.mcl_fold_thres = 11

                    # # /*** mcl_val params ***/
                    if i_metric < 8:
                        self.mcl_inflat = 5.8
                        self.mcl_exp = 6
                    if i_metric >= 8 and i_metric <= 11:
                        self.mcl_inflat = 2.3
                        self.mcl_exp = 4
                    if i_metric > 11:
                        self.mcl_inflat = 1.1
                        self.mcl_exp = 4
                    self.mcl_fold_thres = 11

                    self.umi_len = i_metric
                else:
                    fn_surf = str(i_metric)
                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                fn = self.fn_pref[self.metric] + fn_surf
                if is_trim:
                    self.trim(
                        fastq_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/' + fn,
                        fastq_trimmed_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/trimmed/' + fn,
                        umi_len=self.umi_len,
                    )
                if is_tobam:
                    fas2bam(
                        fastq_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/' + self.comp_cat + '/' + fn + '.fastq.gz',
                        bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/' + self.comp_cat + '/bam/' + fn,
                    ).tobam()
                if is_dedup:
                    # if self.metric == 'seq_errs':
                    #     if i_metric == 0.125 or i_metric == 0.15:
                    #         continue
                    #     else:
                            dedup_ob = dedupPos(
                                mode='internal',
                                method=self.method,
                                # bam_fpn=to('example/data/example.bam'),
                                bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/' + self.comp_cat + '/bam/' + fn + '.bam',
                                pos_tag='PO',
                                mcl_fold_thres=self.mcl_fold_thres,
                                inflat_val=self.mcl_inflat,
                                exp_val=self.mcl_exp,
                                iter_num=100,
                                verbose=False,
                                ed_thres=1,
                                is_sv=False,
                                sv_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/' + self.comp_cat + '/bam/' + fn + '.summary',
                            )
                            dedup_arr.append(dedup_ob.dedup_num)
            df_dedup['pn' + str(i_pn)] = dedup_arr
            print(df_dedup)
        # self.gwriter.generic(
        #     df=df_dedup,
        #     sv_fpn=fastq_fp + self.metric + '/' + str(self.method) + '_' + self.comp_cat + '.txt',
        #     header=True,
        # )

    def trim(self, fastq_fpn, fastq_trimmed_fpn, umi_len):
        trim_params = {
            'read_struct': 'umi_1',
            'umi_1': {
                'len': umi_len,
            },
            'fastq': {
                'fpn': fastq_fpn + '.fastq.gz',
                'trimmed_fpn': fastq_trimmed_fpn + '.fastq.gz',
            },
        }
        umitrim_parser = umitrim(trim_params)
        df = umitrim_parser.todf()
        umitrim_parser.togz(df)
        return 0


if __name__ == "__main__":
    p = umi(
        # metric='pcr_nums',
        # metric='pcr_errs',
        metric='seq_errs',
        # metric='ampl_rates',
        # metric='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        method='directional',
        # method='mcl',
        # method='mcl_val',
        # method='mcl_ed',

        # comp_cat='',
        # comp_cat='ref',
        comp_cat='bipartite',

        # is_trim=True,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,
        fastq_fp=to('data/simu/dimer/pcr8/'),
    )
