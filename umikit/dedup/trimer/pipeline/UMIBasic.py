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
from umikit.dedup.trimer.pipeline import Config
from umikit.dedup.monomer.Relation import relation as umimonorel
from umikit.dedup.monomer.Trace import trace as umimonotrace
from umikit.deduplicate.monomer.DedupPos import dedupPos
from umikit.dedup.monomer.Cluster import cluster as umimonoclust
from umikit.dedup.monomer.Adjacency import adjacency as umitoolmonoadj
from umikit.dedup.monomer.Directional import directional as umitoolmonodirec
from umikit.dedup.monomer.MarkovClustering import markovClustering as umimonomcl
from umikit.dedup.monomer.DBSCAN import dbscan as dbsc
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
                    self.umi_len = self.umi_unit_len_fixed
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
                    self.umi_len = self.umi_unit_len_fixed
                elif self.metric == 'seq_errs':
                    print('=>No.{} sequencing error: {}'.format(id, i_metric))
                    # self.mcl_inflat = 1.1 if i_metric > 0.005 else 2.7
                    # self.mcl_exp = 3
                    self.mcl_fold_thres = 1.6
                    self.mcl_inflat = 1.1 if i_metric > 0.005 else 2.7
                    self.mcl_exp = 2
                    fn_surf = str(id)
                    self.umi_len = self.umi_unit_len_fixed
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

                    self.umi_len = self.umi_unit_len_fixed
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
                    self.umi_len = self.umi_unit_len_fixed
                fn = self.fn_pref[self.metric] + fn_surf
                if is_trim:
                    self.trim(
                        fastq_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/' + fn,
                        fastq_trimmed_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/trimmed/' + fn,
                        umi_len=self.umi_len,
                    )
                if is_tobam:
                    fas2bam(
                        fastq_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/trimmed/' + self.comp_cat + '/' + fn + '.fastq.gz',
                        bam_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/bam/' + self.comp_cat + '/' + fn,
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
                                bam_fpn=fastq_fp + self.metric + '/trimer/permute_' + str(i_pn) + '/bam/' + self.comp_cat + '/' + fn + '.bam',
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
            print(df_dedup)
        self.gwriter.generic(
            df=df_dedup,
            sv_fpn=fastq_fp + self.metric + '/' + str(self.method) + '.txt',
            header=True,
        )

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

    def statistics(self, ):
        return {
            'ccs': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'metric', 'method']),
            },
            'adj': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'metric', 'method']),
            },
            'direc': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'metric', 'method']),
            },
            'mcl_val': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'metric', 'method']),
            },
            'mcl_ed': {
                'df_apv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_disapv_cnt': pd.DataFrame(columns=['0', '1', 'all', 'metric', 'method']),
                'df_apv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_disapv_pct': pd.DataFrame(columns=['0', '1', 'metric', 'method']),
                'df_dedup_cnt': pd.DataFrame(columns=['dedup_cnt', 'metric', 'method']),
            },
        }

    def evaluate(self, ):
        stat = self.statistics()
        for id, i_metric in enumerate(self.metrics[self.metric]):
            if self.metric == 'pcr_nums':
                print('=>at PCR {}'.format(i_metric))
                fastq_fn_surf = str(i_metric)
            elif self.metric == 'pcr_errs':
                print('=>No.{} PCR error: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'seq_errs':
                print('=>No.{} sequencing error: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'ampl_rates':
                print('=>No.{} amplification rate: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'umi_lens':
                print('=>No.{} UMI length: {}'.format(id, i_metric))
                fastq_fn_surf = str(i_metric)
            else:
                fastq_fn_surf = i_metric
            fastq_fn = self.fn_pref[self.metric] + fastq_fn_surf
            umikit = self.umimonorel(
                fastq_path=to('data/simu/') + self.metric + '/trimmed/',
                # fastq_path=to('data/simu/') + self.metric + '/spacer/trimmed/',
                fastq_name=fastq_fn,
                ed_thres=1,
            )

            umitrace = umimonotrace(
                df_fastq=umikit.df_fastq,
                df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                umi_uniq_mapped_rev=umikit.umi_uniq_mapped_rev,
                umi_trace_dict=umikit.umi_trace_dict,
            )

            umitooldirec = umitoolmonodirec()
            umitooladj = umitoolmonoadj()

            ccs_stime = time.time()
            ccs = umimonoclust().cc(graph_adj=umikit.graph_adj)
            ccs_net_num = len([*ccs])
            print('time for obtaining connected components: {time:.3f}s'.format(time=time.time() - ccs_stime))
            print('connected component number: {}'.format(ccs_net_num))

            adj_stime = time.time()
            adj_net_num = umitooladj.umi_tools(
                connected_components=ccs,
                df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                graph_adj=umikit.graph_adj,
            )
            print('time for using Adjacency: {time:.3f}s'.format(time=time.time() - adj_stime))
            print('Adjacency network number: {}'.format(adj_net_num))
            direc_stime = time.time()
            direc_net_num, direc_cc_subs, direc_cc_apvs, direc_cc_disapvs = umitooldirec.umi_tools(
                connected_components=ccs,
                df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                graph_adj=umikit.graph_adj,
            )
            print('time for using Directional: {time:.3f}s'.format(time=time.time() - direc_stime))
            print('Directional network number: {}'.format(direc_net_num))

            # @@@ block: MCL
            mcl_stime = time.time()
            mcl = umimonomcl()
            dbscan = dbsc()
            df_mcl_ccs = mcl.dfclusters(
                connected_components=ccs,
                graph_adj=umikit.graph_adj,
            )
            mcl_num = df_mcl_ccs['clust_num'].sum()
            print('time for using MCL: {time:.3f}s'.format(time=time.time() - mcl_stime))
            print('MCL network number: {}'.format(mcl_num))

            mscmv_val_stime = time.time()
            mscmv_val_len, mscmv_val_clusters, mscmv_val_apv, mscmv_val_disapv = mcl.maxval_val(
                df_mcl_ccs=df_mcl_ccs,
                df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                thres_fold=2,
            )
            mscmv_val_dedup_cnt = mscmv_val_len.sum()
            print('time for using MCL-Val: {time:.3f}s'.format(time=time.time() - mscmv_val_stime))
            print('MCL-Val network number: {}'.format(mscmv_val_dedup_cnt))

            mscmv_ed_stime = time.time()
            mscmv_ed_len, mscmv_ed_clusters, mscmv_ed_apv, mscmv_ed_disapv = mcl.maxval_ed(
                df_mcl_ccs=df_mcl_ccs,
                df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                umi_uniq_mapped_rev=umikit.umi_uniq_mapped_rev,
                thres_fold=1,
            )
            mscmv_ed_dedup_cnt = mscmv_ed_len.sum()
            print('time for using MCL-ED: {time:.3f}s'.format(time=time.time() - mscmv_ed_stime))
            print('MCL-ED network number: {}'.format(mscmv_ed_dedup_cnt))

            dbscan_stime = time.time()
            t = dbscan.dfclusters(
                connected_components=ccs, df_umi_uniq_val_cnt=umikit.df_umi_uniq_val_cnt,
                umi_uniq_mapped_rev=umikit.umi_uniq_mapped_rev,
            )
            # t = dbscan.dfclusters(umikit.umi_uniq_mapped_rev)
            print('time for using DBSCAN: {time:.3f}s'.format(time=time.time() - dbscan_stime))
            print('DBSCAN cluster number: {}'.format(t))

            # print(stat['ccs']['df_dedup_cnt'])
            #
            # stat['ccs']['df_dedup_cnt'].loc[i_metric] = [ccs_net_num] + [str(i_metric)] + ['ccs']
            # stat['adj']['df_dedup_cnt'].loc[i_metric] = [adj_net_num] + [str(i_metric)] + ['adj']

            # df_direc_apv = umitooldirec.formatApvsDisapv(direc_cc_apvs)
            # df_direc_disapv = umitooldirec.formatApvsDisapv(direc_cc_disapvs)
            # stat['direc']['df_apv_cnt'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=df_direc_apv, sort='cnt')) + [str(i_metric)] + ['direc']
            # stat['direc']['df_disapv_cnt'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=df_direc_disapv, sort='cnt')) + [str(i_metric)] + ['direc']
            # stat['direc']['df_apv_pct'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=df_direc_apv, sort='pct')) + [str(i_metric)] + ['direc']
            # stat['direc']['df_disapv_pct'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=mscmv_ed_disapv, sort='pct')) + [str(i_metric)] + ['direc']
            # stat['direc']['df_dedup_cnt'].loc[i_metric] = [direc_net_num] + [str(i_metric)] + ['direc']
            #
            # stat['mcl_val']['df_dedup_cnt'].loc[i_metric] = [mscmv_val_dedup_cnt] + [str(i_metric)] + ['mcl_val']
            #
            # stat['mcl_ed']['df_apv_cnt'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=mscmv_ed_apv, sort='cnt')) + [str(i_metric)] + ['mcl_ed']
            # stat['mcl_ed']['df_disapv_cnt'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=mscmv_ed_disapv, sort='cnt')) + [str(i_metric)] + ['mcl_ed']
            # stat['mcl_ed']['df_apv_pct'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=mscmv_ed_apv, sort='pct')) + [str(i_metric)] + ['mcl_ed']
            # stat['mcl_ed']['df_disapv_pct'].loc[i_metric] = list(umitrace.edgecls(df_list_2d=mscmv_ed_disapv, sort='pct')) + [str(i_metric)] + ['mcl_ed']
            # stat['mcl_ed']['df_dedup_cnt'].loc[i_metric] = [mscmv_ed_dedup_cnt] + [str(i_metric)] + ['mcl_ed']

            # if id == len(self.metrics[self.metric])-1:
            # # if id == 5:
            #     df_dedup_cnt = pd.concat([
            #         # stat['ccs']['df_dedup_cnt'],
            #         # stat['adj']['df_dedup_cnt'],
            #         stat['direc']['df_dedup_cnt'],
            #         stat['mcl_val']['df_dedup_cnt'],
            #         stat['mcl_ed']['df_dedup_cnt'],
            #     ]).reset_index(drop=True)
            #     # print(df_dedup_cnt)
            #     # self.plotv.n1(df_disapv=stat['direc']['df_disapv_cnt'], df_apv=stat['direc']['df_apv_cnt'])
            #     # self.plotv.n1(df_disapv=stat['mcl_ed']['df_disapv_cnt'], df_apv=stat['mcl_ed']['df_apv_cnt'])
            #     self.plotv.n2(df=df_dedup_cnt)
            #     self.plotv.n2dist(df=df_dedup_cnt)
            #     tt[2] = tt1[0]
            #     sns.set_theme(style="ticks")
            #     tt = tt.reset_index(drop=True)
            #     # print(tt.)
            #     # sns.boxplot(x=0, y=1, data=tt,
            #     #             whis=[0, 100], width=.6, palette="vlag")
            #     # sns.stripplot(x=0, y=1, data=tt,
            #     #               size=4, color=".3", linewidth=0)
            #     # sns.jointplot(
            #     #     data=tt,
            #     #     x=0, y=2, hue=1,
            #     #     kind="kde",
            #     #     # palette='Paired'
            #     # )
            #     sd = tt.dropna()
            #
            #     print(sd)
            #     g = sns.lmplot(
            #         data=sd,
            #         x=0, y=2, hue=1,
            #         height=5
            #     )


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
        # method='directional',
        # method='mcl',
        # method='mcl_val',
        method='mcl_ed',

        # comp_cat='ref',
        comp_cat='bipartite',

        # is_trim=True,
        # is_tobam=True,
        # is_dedup=False,

        is_trim=False,
        is_tobam=False,
        is_dedup=True,
        fastq_fp=to('data/simu/umi/'),
    )
    # print(p.evaluate())