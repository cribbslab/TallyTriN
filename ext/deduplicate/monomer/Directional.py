__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"


class directional(object):

    def __init__(self, ):
        pass

    def umi_tools(self, connected_components, df_umi_uniq_val_cnt, graph_adj):
        cc_sub_cnt = []
        cc_subs = {}
        cc_apvs = {}
        cc_disapvs = {}
        for i, cc in connected_components.items():
            cc_sub, apv_node_nbr, disapv_node_nbr = self.umi_tools_(
                df_umi_uniq_val_cnt=df_umi_uniq_val_cnt,
                cc=cc,
                graph_adj=graph_adj,
            )
            # print(cc_sub)
            cc_sub_cnt.append(len(cc_sub))
            cc_subs['cc_' + str(i)] = cc_sub
            cc_apvs['cc_' + str(i)] = apv_node_nbr
            cc_disapvs['cc_' + str(i)] = disapv_node_nbr
        # print(sum(cc_sub_cnt))
        # print(cc_subs)
        # print(cc_apvs)
        # print(cc_disapvs)
        return sum(cc_sub_cnt), cc_subs, cc_apvs, cc_disapvs

    def umi_tools_(self, df_umi_uniq_val_cnt, cc, graph_adj):
        cc_node_sorted = df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()
        nodes = [*cc_node_sorted.keys()]
        # print(nodes)
        node_cp = nodes.copy()
        node_set_remaining = set(node_cp)
        cc_sub = {}
        apv_node_nbr = {}
        disapv_node_nbr = {}
        while nodes:
            e = nodes.pop(0)
            if e in node_set_remaining:
                seen, apv, disapv = self.dfs(e, cc_node_sorted, node_set_remaining, graph_adj)
                cc_sub['node_' + str(e)] = list(seen)
                apv_node_nbr['node_' + str(e)] = apv
                disapv_node_nbr['node_' + str(e)] = disapv
                node_set_remaining = node_set_remaining - seen
                # print('remaining: {}'.format(node_set_remaining))
                # print('disapproval {}'.format(disapv))
            else:
                continue
        # print(disapv_node_nbr)
        return cc_sub, apv_node_nbr, disapv_node_nbr

    def dfs(self, node, node_val_sorted, node_set_remaining, graph_adj):
        visited = set()
        approval = []
        disapproval = []
        g = graph_adj
        def search(node):
            visited.add(node)
            # print(visited)
            for neighbor in g[node]:
                # print(neighbor)
                if neighbor not in visited:
                    if neighbor in node_set_remaining:
                        if node_val_sorted[node] >= 2 * node_val_sorted[neighbor] - 1:
                            approval.append([node, neighbor])
                            search(neighbor)
                        else:
                            disapproval.append([node, neighbor])
        search(node)
        return visited, approval, disapproval

    def formatApvsDisapv(self, cc_dict):
        """
        input format for the Directional method in umi-tools
        :return:
        """
        cc_lvl_df = pd.Series(cc_dict)
        return cc_lvl_df.apply(lambda x: self.dictTo2d(x))

    def dictTo2d(self, x):
        node_list_3d = [*x.values()]
        res_2d = []
        for i in node_list_3d:
            res_2d = res_2d + i
        return res_2d

    def formatCCS(self, cc_dict):
        cc_lvl_df = pd.Series(cc_dict)
        return cc_lvl_df.apply(lambda x: [*x.values()])


if __name__ == "__main__":
    import pandas as pd
    from ext.deduplicate.monomer.Cluster import cluster as umimonoclust

    p = directional()

    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }

    node_val_sorted = pd.Series({
        'A': 456,
        'E': 90,
        'D': 72,
        'B': 2,
        'C': 2,
        'F': 1,
    })
    print(node_val_sorted)
    ccs = umimonoclust().cc(graph_adj=graph_adj)
    print(p.umi_tools(
        connected_components=ccs,
        df_umi_uniq_val_cnt=node_val_sorted,
        graph_adj=graph_adj
    ))

