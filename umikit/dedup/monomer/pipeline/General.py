__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import pandas as pd
from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from simreadflow.read.umi.Filter import filter
from umikit.dedup.monomer.Relation import relation as umimonorel
from umikit.dedup.monomer.Trace import trace as umimonotrace
from umikit.dedup.monomer.Cluster import cluster as umimonoclust
from umikit.dedup.monomer.Adjacency import adjacency as umitoolmonoadj
from umikit.dedup.monomer.Directional import directional as umitoolmonodirec
from umikit.dedup.monomer.MarkovClustering import markovClustering as umimonomcl
from umikit.dedup.monomer.DBSCAN import dbscan as dbsc
from umikit.dedup.monomer.plot.Valid import valid as plotv
from simreadflow.simulate.dispatcher.batch.General import general as simugeneralstarter
from simreadflow.read.similarity.distance.Hamming import hamming
from simreadflow.util.sequence.fastq.Read import read as rfastq
from simreadflow.util.sequence.fastq.Write import write as wfastq
from Path import to


class batch(simugeneralstarter):

    def __init__(self, metric):
        super(batch, self).__init__()
        self.filter = filter()
        self.rfastq = rfastq
        self.wfastq = wfastq
        self.gbfscc = gbfscc()
        self.hamming = hamming()
        self.umimonorel = umimonorel
        self.plotv = plotv()
        # print(self.pcr_errs, self.seq_errs)

        self.metric = metric

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
                print('->at PCR {}'.format(i_metric))
                fastq_fn_surf = str(i_metric)
            elif self.metric == 'pcr_errs':
                print('->No.{} PCR error: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'seq_errs':
                print('->No.{} sequencing error: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'ampl_rates':
                print('->No.{} amplification rate: {}'.format(id, i_metric))
                fastq_fn_surf = str(id)
            elif self.metric == 'umi_lens':
                print('->No.{} UMI length: {}'.format(id, i_metric))
                fastq_fn_surf = str(i_metric)
            else:
                fastq_fn_surf = i_metric
            fastq_fn = self.fastq_fn_pref[self.metric] + fastq_fn_surf
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
    p = batch(
        # metric='pcr_nums',
        # metric='pcr_errs',
        metric='seq_errs',
        # metric='ampl_rates',
        # metric='umi_lens',
    )
    print(p.evaluate())

