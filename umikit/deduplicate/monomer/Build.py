__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import pandas as pd
from Path import to
from umikit.util.Number import number as rannum
from umikit.util.Hamming import hamming
from umikit.network.Edge import edge as guuedge
from umikit.util.Console import console


class build():

    def __init__(self, df, ed_thres, verbose=False):
        """

        Parameters
        ----------
        df
        ed_thres
        """
        self.df = df
        # print(df)
        self.hamming = hamming()
        self.guuedge = guuedge()
        self.console = console()
        self.console.verbose = verbose

        umi_keymap_stime = time.time()
        self.uniq_umis = self.df['umi'].unique()
        self.uniq_umi_num = self.uniq_umis.shape[0]
        self.console.print('===>unique UMI number: {}'.format(self.uniq_umi_num))
        self.umi_uniq_mapped = {k: id for id, k in enumerate(self.uniq_umis)}
        # print(self.umi_uniq_mapped)
        self.umi_uniq_mapped_rev = {id: k for k, id in self.umi_uniq_mapped.items()}
        # print(self.umi_uniq_mapped_rev)
        self.df_umi_uniq_val_cnt = self.df['umi'].value_counts(ascending=False)
        # print(self.df_umi_uniq_val_cnt)
        df_umi_uniq_val_cnt_indexes = self.df_umi_uniq_val_cnt.index
        self.df_umi_uniq_val_cnt.index = [self.umi_uniq_mapped[i] for i in df_umi_uniq_val_cnt_indexes]
        # print(self.df_umi_uniq_val_cnt)
        self.console.print('===>umi keymap time: {:.3f}s'.format(time.time() - umi_keymap_stime))
        self.umi_bam_ids = {}
        for k, v in self.umi_uniq_mapped_rev.items():
            self.umi_bam_ids[k] = df.loc[df['umi'].isin([v])]['id'].values[0]

        ed_list_stime = time.time()
        self.ed_list = self.ed_list_()
        self.df_eds = pd.DataFrame(self.ed_list, columns=['node_1', 'node_2', 'ed'])
        self.df_ed_sel = self.df_eds.loc[self.df_eds['ed'] == ed_thres]
        self.console.print('===>edit distance list construction time: {:.3f}s'.format(time.time() - ed_list_stime))

        edge_list_stime = time.time()
        self.edge_list = self.guuedge.fromdf(self.df_ed_sel, to_tuple=False)
        self.console.print('===>edge list construction time: {:.3f}s'.format(time.time() - edge_list_stime))

        graph_adj_stime = time.time()
        self.graph_adj = {i: [] for i in [*self.umi_uniq_mapped.values()]}
        self.guuedge.graph = self.edge_list
        self.graph_adj_edges = self.guuedge.toAdjacencyDict()
        for k, v in self.graph_adj_edges.items():
            self.graph_adj[k] += v
        # from umikit.deduplicate.monomer.Cluster import cluster as umimonoclust
        # cc = umimonoclust().cc(self.graph_adj)
        # fff = False
        # if len(self.umi_uniq_mapped_rev) > 200:
        #     print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #     print(self.edge_list)
        #     print(self.graph_adj)
        #     print(self.umi_uniq_mapped_rev)
        #     print(self.df_umi_uniq_val_cnt)
        #     fff = True
        # else:
        #     fff = False
            # break
        self.console.print('===>graph adjacency list construction time: {:.3f}s'.format(time.time() - graph_adj_stime))
        self.data_summary = {
            'graph_adj': self.graph_adj,
            'umi_uniq_mapped_rev': self.umi_uniq_mapped_rev,
            'df_umi_uniq_val_cnt': self.df_umi_uniq_val_cnt,
            'umi_bam_ids': self.umi_bam_ids,
            # 'fff': fff,
        }

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
        """

        Parameters
        ----------
        x

        Returns
        -------

        """
        c = x.split('_')[0].split('-')
        if c[1] == 'init':
            return -1
        else:
            return c[2]




if __name__ == "__main__":
    umikit = build(
        # bam_fpn=to('example/data/example.bam'),
        bam_fpn=to('example/data/RM82CLK1_S1_XT_gene.bam'),
        ed_thres=5,
    )