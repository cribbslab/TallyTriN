__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from src.util.sequence.fastq.Read import read as rfastq
from src.util.sequence.fastq.Write import write as wfastq
from src.sequencing.reads.similarity.distance.Hamming import hamming
from src.util.graph.undirected.unweighted.repr.Edge import edge as guuedge


class relation(object):

    def __init__(self, fastq_path, fastq_name, ed_thres):
        self.rfastq = rfastq
        self.wfastq = wfastq
        self.hamming = hamming()
        self.guuedge = guuedge()
        print('===>read the fastq file {}.fastq.gz'.format(fastq_name))
        read_stime = time.time()
        self.names, self.seqs, _, _ = self.rfastq().fromgz(
            fastq_path=fastq_path,
            fastq_name=fastq_name,
            method='pyfastx',
        )
        print('===>file read time: {:.3f}s'.format(time.time() - read_stime))

        umi_df_stime = time.time()
        self.df_fastq = pd.DataFrame(self.seqs, columns=['seq_raw'])
        self.df_fastq['name'] = self.names
        self.df_fastq['umi'] = self.df_fastq['name'].apply(lambda x: x.split('_')[1])
        self.df_fastq['umi#'] = self.df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[0])
        self.df_fastq['umi_src'] = self.df_fastq['name'].apply(lambda x: x.split('_')[0].split('-')[1])
        self.df_fastq['umi_pcr#'] = self.df_fastq['name'].apply(lambda x: self.pcrnum(x))
        print('===>umi to df time: {:.3f}s'.format(time.time() - umi_df_stime))
        # print(df['umi#'])

        umi_keymap_stime = time.time()
        self.uniq_umis = self.df_fastq['umi'].unique()
        self.uniq_umi_num = self.uniq_umis.shape[0]
        print('===>unique UMI number: {}'.format(self.uniq_umi_num))
        self.umi_uniq_mapped = {k: id for id, k in enumerate(self.uniq_umis)}
        self.umi_uniq_mapped_rev = {id: k for k, id in self.umi_uniq_mapped.items()}
        self.df_umi_uniq_val_cnt = self.df_fastq['umi'].value_counts()
        # print(self.df_umi_uniq_val_cnt)
        df_umi_uniq_val_cnt_indexes = self.df_umi_uniq_val_cnt.index
        self.df_umi_uniq_val_cnt.index = [self.umi_uniq_mapped[i] for i in df_umi_uniq_val_cnt_indexes]
        print('===>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))

        umi_trace_dict_stime = time.time()
        self.umi_trace_dict = {}
        self.df_fastq_umi_gp = self.df_fastq.groupby(['umi'])
        for uniq_umi in self.uniq_umis:
            self.umi_trace_dict[self.umi_uniq_mapped[uniq_umi]] = self.df_fastq_umi_gp.get_group(uniq_umi)['umi#'].unique().astype(np.int).tolist()
        # print(self.umi_trace)
        print('===>umi trace dict time: {:.3f}s'.format(time.time() - umi_trace_dict_stime))

        ed_list_stime = time.time()
        self.ed_list = self.ed_list_()
        self.df_eds = pd.DataFrame(self.ed_list, columns=['node_1', 'node_2', 'ed'])
        self.df_ed_sel = self.df_eds.loc[self.df_eds['ed'] == ed_thres]
        print('===>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))

        edge_list_stime = time.time()
        self.edge_list = self.guuedge.fromdf(self.df_ed_sel, to_tuple=False)
        print('===>edge list construction time: {:.3f}s'.format(time.time() - edge_list_stime))

        graph_adj_stime = time.time()
        self.graph_adj = {i: [] for i in [*self.umi_uniq_mapped.values()]}
        self.guuedge.graph = self.edge_list
        self.graph_adj_edges = self.guuedge.toAdjacencyDict()
        for k, v in self.graph_adj_edges.items():
            self.graph_adj[k] += v
        print('===>graph adjacency construction time: {:.3f}s'.format(time.time() - graph_adj_stime))

    def ed_list_(self, ):
        eds = []
        for i in range(self.uniq_umi_num):
            for j in range(i + 1, self.uniq_umi_num):
                l = self.uniq_umis[i]
                r = self.uniq_umis[j]
                # if self.umi_uniq_mapped[l] == 31:
                #     print(l)
                # if self.umi_uniq_mapped[r] == 50:
                #     print(r)
                eds.append([
                    self.umi_uniq_mapped[l],
                    self.umi_uniq_mapped[r],
                    self.hamming.general(l, r),
                ])
        # print(len(eds))
        return eds

    def pcrnum(self, x):
        c = x.split('_')[0].split('-')
        if c[1] == 'init':
            return -1
        else:
            return c[2]