import time
import numpy as np


class trace(object):

    def __init__(self, df_fastq, df_umi_uniq_val_cnt, umi_uniq_mapped_rev, umi_trace_dict):
        self.df_fastq = df_fastq
        self.umi_uniq_mapped_rev = umi_uniq_mapped_rev
        self.df_umi_uniq_val_cnt = df_umi_uniq_val_cnt
        self.umi_trace_dict = umi_trace_dict

    def edgecls(self, df_list_2d, sort='cnt'):
        stime = time.time()
        print('------>edges to be traced is {}'.format(df_list_2d.shape[0]))
        res = df_list_2d.apply(
            lambda list_2d: self.by01(list_2d),
        )
        # print(res)
        if sort == 'cnt':
            total = res.apply(lambda x: len(x))
            total_0 = res.apply(lambda x: sum(1 for ele in x if ele == 0))
            total_1 = res.apply(lambda x: sum(1 for ele in x if ele == 1))
            total_lens = total.sum()
            total_0_lens = total_0.sum()
            total_1_lens = total_1.sum()
            print('------>trace edge cls time {time:.2f}s'.format(time=time.time()-stime))
            return total_0_lens, total_1_lens, total_lens
        elif sort == 'pct':
            total_0_pcts = res.apply(lambda x: sum(1 for ele in x if ele == 0) / len(x) if len(x) != 0 else None)
            total_1_pcts = res.apply(lambda x: sum(1 for ele in x if ele == 1) / len(x) if len(x) != 0 else None)
            # print(total_0_pcts.loc[total_0_pcts != None])
            # print(total_1_pcts.loc[total_1_pcts != None])
            return total_0_pcts, total_1_pcts

    def by01(self, list_2d):
        """

        :param list_2d: a 2d list, one with 2-sized vector consisting of 2 nodes.
        :return:
        """
        trace_marks = []
        for nodes in list_2d:
            umi_ori_node_1 = self.umi_trace_dict[nodes[0]]
            umi_ori_node_2 = self.umi_trace_dict[nodes[1]]
            # if len(umi_ori_node_1) != 1 and len(umi_ori_node_2) != 1:
            #     trace_marks.append(0)
            # else:
            intxn = set(umi_ori_node_1) & set(umi_ori_node_2)
            # print('l {} r {}, intersection {}'.format(umi_ori_node_1, umi_ori_node_2, intxn))
            if len(intxn) != 0:
                trace_marks.append(1)
            else:
                trace_marks.append(0)
        return trace_marks

    def matchRepresentative(self, df):
        print('------>representative to be traced is {}'.format(df.shape[0]))
        stime = time.time()
        res = df.apply(
            lambda list_2d: self.maxval(list_2d),
        )
        # print(res)
        ttt = set()
        total1 = res.values
        for i in total1:
            ttt.update(i)
            # ttt = ttt + i
        total = res.apply(lambda x: len(x))
        # print(total)
        total_len = total.sum()
        # print(total_len)
        print('------>trace representative time {time:.2f}s'.format(time=time.time() - stime))
        return len(ttt)/50

    def maxval(self, list_2d):
        repr_max_nodes = []
        # repr_max_nodes = set()
        for sub_cc_nodes in list_2d:
            node_val_sorted = self.df_umi_uniq_val_cnt.loc[
                self.df_umi_uniq_val_cnt.index.isin(sub_cc_nodes)
            ].sort_values(ascending=False).to_dict()
            max_node = [*node_val_sorted.keys()][0]
            umi_ori_max_nodes = self.umi_trace_dict[max_node]
            if len(umi_ori_max_nodes) != 1:
                repr_max_nodes = repr_max_nodes + []
            else:
                repr_max_nodes = repr_max_nodes + umi_ori_max_nodes
        # print(len(repr_max_nodes))
        return repr_max_nodes