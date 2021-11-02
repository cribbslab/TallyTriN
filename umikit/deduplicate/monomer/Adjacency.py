__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"


class adjacency(object):

    def umi_tools(self, connected_components, df_umi_uniq_val_cnt, graph_adj):
        """
        Examples
        --------
        umi_tools adjacency wrap

        Parameters
        ----------
        connected_components
        df_umi_uniq_val_cnt
        graph_adj

        Returns
        -------

        """
        tcl = []
        cc_subs = {}
        for i, cc in connected_components.items():
            # print('cc: ', cc)
            step, cc_sub = self.umi_tools_(df_umi_uniq_val_cnt=df_umi_uniq_val_cnt, cc=cc, graph_adj=graph_adj)
            tcl.append(step)
            cc_subs['cc_' + str(i)] = cc_sub
            # print(self.umi_tools(cc_sorted))
        return {
            'count': sum(tcl),
            'clusters': cc_subs,
        }

    def umi_tools_(self, df_umi_uniq_val_cnt, cc, graph_adj):
        """
        umi_tools adjacency

        Parameters
        ----------
        df_umi_uniq_val_cnt
            unique umi counts
        cc
            connected_components
        graph_adj
            the adjacency list of a graph

        Returns
        -------

        """
        cc_umi_sorted = df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict()
        cc_sorted = [*cc_umi_sorted.keys()]
        visited = set()
        step = 1
        subcomponents = {}
        cc_set = set(cc_sorted)
        while cc_sorted:
            e = cc_sorted.pop(0)
            subcomponents[e] = []
            for node in graph_adj[e]:
                if node not in visited:
                    subcomponents[e].append(node)
            visited.add(e)
            visited.update(graph_adj[e])
            subcomponents[e] = graph_adj[e]
            # print('the ccurent ele popping out: {} {}'.format(e, visited))
            if visited == cc_set:
                # print(step)
                break
            else:
                step += 1
        vertex = [*subcomponents.keys()]
        cc_sub = {}
        for k, v in subcomponents.items():
            cc_sub['node_' + str(k)] = [k]
            for i in v:
                if i not in vertex:
                    cc_sub['node_' + str(k)].append(i)
        return step, cc_sub

    def decompose(self, cc_sub_dict):
        """

        Parameters
        ----------
        cc_sub_dict

        Returns
        -------

        """
        cc_cnt = 0
        ccs = {}
        for k1, v1 in cc_sub_dict.items():
            for k2, v2 in v1.items():
                ccs[cc_cnt] = v2
                cc_cnt += 1
        return ccs


if __name__ == "__main__":
    import pandas as pd
    from umikit.deduplicate.monomer.Cluster import cluster as umimonoclust

    p = adjacency()

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

