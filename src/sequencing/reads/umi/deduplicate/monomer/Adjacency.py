__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"


class adjacency(object):

    def umi_tools(self, connected_components, df_umi_uniq_val_cnt, graph_adj):
        tcl = []
        for i, cc in connected_components.items():
            # print('cc: ', cc)
            tcl.append(self.umi_tools_(df_umi_uniq_val_cnt=df_umi_uniq_val_cnt, cc=cc, graph_adj=graph_adj))
            # print(self.umi_tools(cc_sorted))
        return sum(tcl)

    def umi_tools_(self, df_umi_uniq_val_cnt, cc, graph_adj):
        cc_umi_sorted = df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()
        cc_sorted = [*cc_umi_sorted.keys()]
        visited = set()
        step = 1
        cc_set = set(cc_sorted)
        while cc_sorted:
            e = cc_sorted.pop(0)
            visited.add(e)
            visited.update(graph_adj[e])
            # print('the ccurent ele popping out: {} {}'.format(e, visited))
            if visited == cc_set:
                # print(step)
                break
            else:
                step += 1
        return step


if __name__ == "__main__":
    p = adjacency()