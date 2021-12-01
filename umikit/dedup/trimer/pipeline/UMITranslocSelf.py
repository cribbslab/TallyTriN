__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from umikit.fastq.Convert import convert as fas2bam
from umikit.trim.Template import template as umitrim
from umikit.util.Writer import writer as gwriter
from simreadflow.util.file.read.Reader import reader as gfreader
from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from umikit.dedup.trimer.pipeline import Config
from umikit.dedup.monomer.Relation import relation as umimonorel
from umikit.deduplicate.monomer.DedupPos import dedupPos
# from umikit.plot.Valid import valid as plotv
from simreadflow.util.random.Number import number as rannum
from Path import to
from umikit.align.Read import read as aliread
import textwrap
from simreadflow.read.similarity.distance.Hamming import hamming


class umiTranslocSelf(Config.config):

    def __init__(self, metric, method, comp_cat, umi_lib_fp=None, fastq_fp=None, is_trim=False, is_tobam=False, is_self_healing=False, is_dedup=False):
        super(umiTranslocSelf, self).__init__()
        self.metric = metric
        self.comp_cat = comp_cat
        self.gbfscc = gbfscc()
        self.gfreader = gfreader()
        self.rannum = rannum()
        self.umimonorel = umimonorel
        self.umi_lib_fp = umi_lib_fp
        self.method = method
        self.gwriter = gwriter()
        self.seq_num = 100
        df_dedup = pd.DataFrame()
        df_proc_bam = pd.DataFrame()

        for i_pn in range(self.permutation_num):
            self.df_umi = self.gfreader.generic(df_fpn=self.umi_lib_fp + self.metric + '/permute_' + str(i_pn) + '/umi.txt')[0].values
            self.df_umi = pd.DataFrame.from_dict({i: e for i, e in enumerate(self.df_umi)}, orient='index', columns=['raw'])
            self.df_umi['index'] = self.df_umi.index
            # print(self.df_umi[0].tolist())
            self.df_umi['collap'] = self.df_umi['raw'].apply(lambda x: ''.join([i[0] for i in textwrap.wrap(x, 3)]))
            self.umi_collap_map = {i: e for i, e in enumerate(self.df_umi['collap'])}
            # print(self.umi_collap_map)
            # reads = np.reshape(self.df_umi[['collap', 'index']].values.tolist(), (int(self.seq_num / 2), 4))
            # print(pd.DataFrame(reads))

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
                    self.mcl_fold_thres = 1.6
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
                        umi_len=self.umi_len * 2,
                    )
                if is_tobam:
                    fas2bam(
                        fastq_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/trimmed/' + fn + '.fastq.gz',
                        bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/bam/' + fn,
                    ).tobam()
                if is_self_healing:
                    self.alireader = aliread(
                        bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/bam/' + fn + '.bam',
                        verbose=False
                    )
                    self.df_bam = self.alireader.todf(tags=['PO'])
                    self.df_bam['names'] = self.df_bam.query_name.apply(lambda x: self.bamproc(x))
                    self.df_bam['r1_id'] = self.df_bam.names.apply(lambda x: x[0])
                    self.df_bam['r2_id'] = self.df_bam.names.apply(lambda x: x[1])
                    self.df_bam['transloc_stat'] = self.df_bam.names.apply(lambda x: x[2])
                    self.df_bam['transloc_side'] = self.df_bam.names.apply(lambda x: x[3])
                    self.df_bam['sam_id'] = self.df_bam.names.apply(lambda x: x[4])
                    self.df_bam['source'] = self.df_bam.names.apply(lambda x: x[5])
                    self.df_bam['umi_l'] = self.df_bam.names.apply(lambda x: x[6][:36])
                    self.df_bam['umi_r'] = self.df_bam.names.apply(lambda x: x[6][36:])
                    self.df_bam['umi_corr_r'] = self.df_bam.umi_r.apply(lambda x: self.correct(x))
                    self.df_fake = self.df_bam.loc[(self.df_bam['transloc_stat'] == 'fake_yes')]
                    # self.df_fake['r2_umi_ref'] = self.df_fake.r1_id.apply(lambda x: print(int(x) + 1))
                    self.df_fake['r2_umi_ref'] = self.df_fake.r1_id.apply(lambda x: self.umi_collap_map[int(x) + 1])
                    # print(self.df_fake['r2_umi_ref'])
                    self.df_fake['cons'] = self.df_fake.apply(lambda x: self.screen1(x), axis=1)
                    self.df_fake['hamm'] = self.df_fake.apply(lambda x: self.screen2(x), axis=1)
                    # print(self.df_fake['hamm'])
                    # print(self.df_fake.loc[self.df_fake['cons'] == 0][['umi_corr_r', 'r2_umi_ref', 'r1_id', 'r2_id', 'transloc_stat', 'source', 'sam_id']])
                    print(self.df_fake.loc[self.df_fake['hamm'] > 6].shape)
                if is_dedup:
                    # if self.metric == 'seq_errs':
                    #     if i_metric == 0.125 or i_metric == 0.15:
                    #         continue
                    #     else:
                            dedup_ob = dedupPos(
                                mode='internal',
                                method=self.method,
                                # bam_fpn=to('example/data/example.bam'),
                                bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/bam/' + fn + '.bam',
                                pos_tag='PO',
                                mcl_fold_thres=self.mcl_fold_thres,
                                inflat_val=self.mcl_inflat,
                                exp_val=self.mcl_exp,
                                iter_num=100,
                                verbose=False,
                                ed_thres=1,
                                is_sv=False,
                                sv_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/summary/' + self.comp_cat + '/' + fn,
                            )
                            dedup_arr.append(dedup_ob.dedup_num)
            df_dedup['pn' + str(i_pn)] = dedup_arr
            # print(df_dedup)
        self.gwriter.generic(
            df=df_dedup,
            sv_fpn=fastq_fp + self.metric + '/' + str(self.method) + '_' + self.comp_cat + '.txt',
            header=True,
        )

    def screen(self, x):
        if x['r1_id'] != x['r2_id']:
            return 1
        else:
            return 0

    def screen1(self, x):
        if x['umi_corr_r'] != x['r2_umi_ref']:
            return 1
        else:
            return 0

    def screen2(self, x):
        return hamming().general(x['umi_corr_r'], x['r2_umi_ref'])

    def bamproc(self, x):
        t = x.split('-')
        # print(t)
        # print(x)
        tt = t[5].split('_')
        # print(t[0], t[1], t[2], t[3], t[4], tt[0], tt[1])
        return t[0], t[1], t[2], t[3], t[4], tt[0], tt[1]

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

    def bamprocSlow(self, x):
        """
                self.df_bam[
                        ['read', 'r1_id', 'r2_id', 'transloc_stat', 'transloc_side', 'sam_id', 'source']
                    ] = self.df_bam.apply(lambda x: self.bamproc(x), axis=1)
        :param x:
        :return:
        """
        t = x['query_name'].split('-')
        # print(t)
        # print(x)
        tt = t[5].split('_')
        # print(t[0], t[1], t[2], t[3], t[4], tt[0], tt[1])
        return pd.Series({
            'read': t[0],
            'r1_id': t[1],
            'r2_id': t[2],
            'transloc_stat': t[3],
            'transloc_side': t[4],
            'sam_id': tt[0],
            'source': tt[1],
        })

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
    p = umiTranslocSelf(
        # metric='pcr_nums',
        # metric='pcr_errs',
        metric='seq_errs',
        # metric='ampl_rates',
        # metric='umi_lens',

        # method='unique',
        # method='cluster',
        # method='adjacency',
        # method='directional',
        # method='mcl',
        method='mcl_val',
        # method='mcl_ed',

        # comp_cat='ref',
        comp_cat='bipartite',

        # is_trim=True,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_self_healing=True,
        is_dedup=False,
        umi_lib_fp=to('data/simu/transloc/trimer/single_read/'),
        fastq_fp=to('data/simu/transloc/trimer/single_read/'),
    )
    # print(p.evaluate())