__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from umikit.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
import networkx as nx


class cluster(object):

    def __init__(self, ):
        pass

    def cc(self, graph_adj):
        connected_components = list(gbfscc().deque(graph_adj))
        return {i: cc for i, cc in enumerate(connected_components)}

    def ccnx(self, edge_list):
        # import matplotlib.pyplot as plt
        G = nx.Graph()
        for edge in edge_list:
            G.add_edge(edge[0], edge[1])
        # nx.draw_networkx(G=G, with_labels=False)
        # plt.show()
        return {i: G.subgraph(cc).nodes() for i, cc in enumerate(nx.connected_components(G))}


if __name__ == "__main__":
    p = cluster()

    graph_adj = {
        'A': ['B', 'C'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }

    graph_umi_tools = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }

    print(p.cc(graph_umi_tools))

