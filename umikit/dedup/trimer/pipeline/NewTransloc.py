__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import numpy as np
import json
import pandas as pd
pd.options.mode.chained_assignment = None
from umikit.util.Writer import writer as gwriter
from simreadflow.util.file.read.Reader import reader as gfreader
from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from umikit.dedup.trimer.pipeline import Config
from umikit.dedup.monomer.Relation import relation as umimonorel
from simreadflow.util.random.Number import number as rannum
from Path import to
from umikit.align.Read import read as aliread
import textwrap
from simreadflow.read.similarity.distance.Hamming import hamming


class newTransloc(Config.config):

    def __init__(self, metric, section, tc_thres=None, umi_lib_fp=None, fastq_fp=None, fn_suffix='', sv_cnt_lib_fpn=''):
        super(newTransloc, self).__init__()
        self.metric = metric
        self.tc_thres = tc_thres
        self.fn_suffix = fn_suffix
        self.gbfscc = gbfscc()
        self.gfreader = gfreader()
        self.gwriter = gwriter()
        self.rannum = rannum()
        self.umimonorel = umimonorel
        self.umi_lib_fp = umi_lib_fp
        self.sv_cnt_lib_fpn = sv_cnt_lib_fpn
        self.section = section
        self.seq_num = 100

        self.real_yes = []
        self.real_no = []
        self.fake_yes = []
        self.chi_yes = []
        self.cnt_dict = {}
        for i_pn in range(self.permutation_num):
            self.df_umi = self.gfreader.generic(df_fpn=self.umi_lib_fp + self.metric + '/permute_' + str(i_pn) + '/umi.txt')[0].values
            self.df_umi = pd.DataFrame.from_dict({i: e for i, e in enumerate(self.df_umi)}, orient='index', columns=['raw'])
            self.df_umi['index'] = self.df_umi.index
            # print(self.df_umi[0].tolist())
            self.df_umi['collap'] = self.df_umi['raw'].apply(lambda x: ''.join([i[0] for i in textwrap.wrap(x, 3)]))
            self.umi_collap_map = {i: e for i, e in enumerate(self.df_umi['collap'])}
            print(self.umi_collap_map)
            # reads = np.reshape(self.df_umi[['collap', 'index']].values.tolist(), (int(self.seq_num / 2), 4))
            # print(pd.DataFrame(reads))
            self.pn_real_yes = []
            self.pn_real_no = []
            self.pn_fake_yes = []
            self.pn_chi_yes = []
            self.cnt_dict[i_pn] = {}
            if self.tc_thres == -1:
                self.tc_thres = self.cnt_paired_umis.quantile(.1)
            for id, i_metric in enumerate(self.metric_vals[self.metric]):
                self.cnt_dict[i_pn][i_metric] = {}

                if self.metric == 'pcr_nums':
                    print('=>at PCR {}'.format(i_metric))
                    fn_surf = str(i_metric)
                    self.mcl_inflat = 2.3
                    self.mcl_exp = 2
                    self.mcl_fold_thres = 1
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
                else:
                    fn_surf = str(i_metric)
                    self.umi_len = self.umi_unit_len_fixed * self.umi_unit_pattern
                fn = self.fn_pref[self.metric] + fn_surf

                self.alireader = aliread(
                    bam_fpn=fastq_fp + self.metric + '/permute_' + str(i_pn) + '/bam/' + fn + '.bam',
                    verbose=False,
                )
                self.df_bam = self.alireader.todf(tags=['PO'])
                self.df_bam['names'] = self.df_bam['query_name'].apply(lambda x: self.bamproc(x))
                self.df_bam['r1_id'] = self.df_bam['names'].apply(lambda x: x[0])
                self.df_bam['r2_id'] = self.df_bam['names'].apply(lambda x: x[1])
                self.df_bam['transloc_stat'] = self.df_bam['names'].apply(lambda x: x[2])
                self.df_bam['transloc_side'] = self.df_bam['names'].apply(lambda x: x[3])
                # print(self.df_bam[self.df_bam['transloc_side'] != 'none'])
                self.df_bam['sam_id'] = self.df_bam['names'].apply(lambda x: x[4])
                self.df_bam['source'] = self.df_bam['names'].apply(lambda x: x[5])
                self.df_bam['umi_l'] = self.df_bam['names'].apply(lambda x: x[6][:30])
                self.df_bam['umi_r'] = self.df_bam['names'].apply(lambda x: x[6][30:])
                self.df_bam['umi'] = self.df_bam.apply(lambda x: x['umi_l'] + x['umi_r'], axis=1)

                self.df_bam['umi_corr_l'] = self.df_bam['umi_l'].apply(lambda x: self.correct(x))
                self.df_bam['umi_corr_r'] = self.df_bam['umi_r'].apply(lambda x: self.correct(x))
                self.df_bam['umi_corr'] = self.df_bam.apply(lambda x: x['umi_corr_l'] + x['umi_corr_r'], axis=1)
                # print(self.df_bam['umi_corr'])
                # print(self.df_bam['umi_corr_l'])
                # print(self.df_bam['umi_corr_r'])

                if self.section == 'dc_by_cnt':
                    self.cnt_paired_umis = self.df_bam['umi'].value_counts()
                    self.hash_paired_umis = self.cnt_paired_umis.to_dict()
                    self.df_bam['chimeric_mark'] = self.df_bam['umi_corr'].apply(
                        lambda x: 1 if self.hash_paired_umis[x] > self.tc_thres else 0
                    )
                    self.chimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 0]
                    self.nonchimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 1]
                    fake_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'fake_yes'].shape[0]
                    real_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_yes'].shape[0]
                    real_no_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_no'].shape[0]
                    chi_num = len(self.chimerical_ids['chimeric_mark'].index)
                    self.pn_fake_yes.append(fake_yes_num)
                    self.pn_real_yes.append(real_yes_num)
                    self.pn_real_no.append(real_no_num)
                    self.pn_chi_yes.append(chi_num)
                    print(fake_yes_num)
                    print(real_yes_num)
                    print(real_no_num)
                    print(chi_num)

                if self.section == 'dc_by_vote':
                    self.cnt_left_umis = self.df_bam['umi_l'].value_counts()
                    self.cnt_right_umis = self.df_bam['umi_r'].value_counts()
                    self.cnt_paired_umis = self.df_bam['umi'].value_counts()
                    self.hash_left_umis = self.cnt_left_umis.to_dict()
                    self.hash_right_umis = self.cnt_right_umis.to_dict()
                    self.hash_paired_umis = self.cnt_paired_umis.to_dict()

                    self.df_bam['chimeric_mark'] = self.df_bam['umi'].apply(lambda x: self.vote(
                        umi=x,
                        hash_left_umis=self.hash_left_umis,
                        hash_right_umis=self.hash_right_umis,
                        hash_paired_umis=self.hash_paired_umis,
                    ))

                    self.chimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 0]
                    self.nonchimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 1]

                    fake_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'fake_yes'].shape[0]
                    real_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_yes'].shape[0]
                    real_no_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_no'].shape[0]
                    chi_num = len(self.chimerical_ids['chimeric_mark'].index)
                    self.pn_fake_yes.append(fake_yes_num)
                    self.pn_real_yes.append(real_yes_num)
                    self.pn_real_no.append(real_no_num)
                    self.pn_chi_yes.append(chi_num)

                if self.section == 'dc_control':
                    self.cnt_l_umis = self.df_bam['umi_l'].value_counts()
                    self.hash_paired_umis = self.cnt_l_umis.to_dict()
                    self.df_bam['chimeric_mark'] = self.df_bam['umi_l'].apply(
                        lambda x: 1 if self.hash_paired_umis[x] > self.tc_thres else 0
                    )
                    self.chimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 0]
                    self.nonchimerical_ids = self.df_bam.loc[self.df_bam['chimeric_mark'] == 1]
                    fake_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'fake_yes'].shape[0]
                    real_yes_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_yes'].shape[0]
                    real_no_num = self.chimerical_ids[self.chimerical_ids['transloc_stat'] == 'real_no'].shape[0]
                    chi_num = len(self.chimerical_ids['chimeric_mark'].index)
                    self.pn_fake_yes.append(fake_yes_num)
                    self.pn_real_yes.append(real_yes_num)
                    self.pn_real_no.append(real_no_num)
                    self.pn_chi_yes.append(chi_num)

                if self.section == 'cnt_lib':
                    self.l_cnt = self.df_bam['umi_l'].value_counts()
                    self.r_cnt = self.df_bam['umi_r'].value_counts()
                    self.lr_cnt = self.df_bam['umi'].value_counts()
                    self.l_corr_cnt = self.df_bam['umi_corr_l'].value_counts()
                    self.r_corr_cnt = self.df_bam['umi_corr_r'].value_counts()
                    self.lr_corr_cnt = self.df_bam['umi_corr'].value_counts()
                    # print(self.l_cnt.values)
                    # print(self.r_cnt.values)
                    # print(self.lr_cnt.values)
                    # print(self.l_corr_cnt.values)
                    # print(self.r_corr_cnt.values)
                    # print(self.lr_corr_cnt.values)
                    self.cnt_dict[i_pn][i_metric]['l_cnt'] = self.l_cnt.values.tolist()
                    self.cnt_dict[i_pn][i_metric]['r_cnt'] = self.r_cnt.values.tolist()
                    self.cnt_dict[i_pn][i_metric]['lr_cnt'] = self.lr_cnt.values.tolist()
                    self.cnt_dict[i_pn][i_metric]['l_corr_cnt'] = self.l_corr_cnt.values.tolist()
                    self.cnt_dict[i_pn][i_metric]['r_corr_cnt'] = self.r_corr_cnt.values.tolist()
                    self.cnt_dict[i_pn][i_metric]['lr_corr_cnt'] = self.lr_corr_cnt.values.tolist()
            # print(self.cnt_dict)
            with open(self.sv_cnt_lib_fpn, 'w') as fp:
                json.dump(self.cnt_dict, fp)
            with open(self.sv_cnt_lib_fpn) as fp1:
                self.ft_n30_dict = json.load(fp1)
            print(self.ft_n30_dict)
            self.real_yes.append(self.pn_real_yes)
            self.real_no.append(self.pn_real_no)
            self.fake_yes.append(self.pn_fake_yes)
            self.chi_yes.append(self.pn_chi_yes)

        # self.gwriter.generic(
        #     df=self.real_yes,
        #     sv_fpn=fastq_fp + self.metric + '/real_yes_' + self.section + '_' + str(tc_thres) + '_' + self.fn_suffix + '.txt',
        #     header=True,
        # )
        # self.gwriter.generic(
        #     df=self.real_no,
        #     sv_fpn=fastq_fp + self.metric + '/real_no_' + self.section + '_' + str(tc_thres) + '_' + self.fn_suffix + '.txt',
        #     header=True,
        # )
        # self.gwriter.generic(
        #     df=self.fake_yes,
        #     sv_fpn=fastq_fp + self.metric + '/fake_yes_' + self.section + '_' + str(tc_thres) + '_' + self.fn_suffix + '.txt',
        #     header=True,
        # )
        # self.gwriter.generic(
        #     df=self.chi_yes,
        #     sv_fpn=fastq_fp + self.metric + '/chi_yes_' + self.section + '_' + str(tc_thres) + '_' + self.fn_suffix + '.txt',
        #     header=True,
        # )

    def vote(self, umi, hash_left_umis, hash_right_umis, hash_paired_umis):
        len_l = len(next(iter(hash_left_umis)))
        len_r = len(next(iter(hash_right_umis)))
        cnt_paired = hash_paired_umis[umi]

        cnt_l = hash_left_umis[umi[:len_l]]
        cnt_r = hash_right_umis[umi[len_l:len_l+len_r]]
        if cnt_paired <= self.tc_thres:
            if cnt_l < cnt_paired or cnt_r > cnt_paired:
                return 1
            else:
                return 0
        else:
            return 1

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
        vernier = [i for i in range(30) if i % 3 == 0]
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


if __name__ == "__main__":
    p = newTransloc(
        metric='seq_errs',
        # section='dc_by_cnt',
        # section='dc_by_vote',
        # section='dc_control',
        section='cnt_lib',
        # section='plot',

        tc_thres=5,
        # fn_suffix='corr',
        fn_suffix='sm&lg',
        umi_lib_fp=to('data/simu/transloc/trimer/single_read/pcr8_umi30/'),
        fastq_fp=to('data/simu/transloc/trimer/single_read/pcr8_umi30/'),
        sv_cnt_lib_fpn=to('data/simu/transloc/trimer/single_read/pcr8_umi30/cnt_lib.json'),
    )