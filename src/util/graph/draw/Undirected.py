import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx


class undirected(object):

    def __init__(self, ):
        pass

    def adaa(self, ):
        import matplotlib.pyplot as plt
        import networkx as nx

        G = nx.Graph()

        G.add_edge('a', 'b', weight=0.6)
        G.add_edge('a', 'c', weight=0.2)
        G.add_edge('c', 'd', weight=0.1)
        G.add_edge('c', 'e', weight=0.7)
        G.add_edge('c', 'f', weight=0.9)
        G.add_edge('a', 'd', weight=0.3)

        elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0.5]
        esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]

        pos = nx.spring_layout(G)  # positions for all nodes

        # nodes
        nx.draw_networkx_nodes(G, pos, node_size=700)

        # edges
        nx.draw_networkx_edges(G, pos, edgelist=elarge,
                               width=6)
        nx.draw_networkx_edges(G, pos, edgelist=esmall,
                               width=6, alpha=0.5, edge_color='b', style='dashed')

        # labels
        nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

        plt.axis('off')
        plt.savefig("weight.jpg")
        plt.show()


    def plot_mute_graph(self, ):
        """Plot MuTE example network.

        Network of 5 AR-processes, which is used as an example the paper
        on the MuTE toolbox (Montalto, PLOS ONE, 2014, eq. 14). The
        network consists of five autoregressive (AR) processes with model
        orders 2 and les and the following (non-linear) couplings:

            >>> 0 -> 1, u = 2
            >>> 0 -> 2, u = 3
            >>> 0 -> 3, u = 2 (non-linear)
            >>> 3 -> 4, u = 1
            >>> 4 -> 3, u = 1

        Returns:
            Figure handle
                Figure object from the matplotlib package
        """
        graph = nx.DiGraph()
        graph.add_nodes_from(np.arange(5))
        # graph.add_edges_from([(0, 1), (0, 2), (0, 3), (3, 4), (4, 3)])
        graph.add_weighted_edges_from([(0, 1, 2), (0, 2, 3), (0, 3, 2), (3, 4, 1),
                                       (4, 3, 1)], weight='delay')
        pos = {
            0: np.array([1, 1]),
            1: np.array([0, 2]),
            2: np.array([0, 0]),
            3: np.array([2, 1]),
            4: np.array([3, 1]),
        }
        fig = plt.figure()
        nx.draw(graph, pos=pos, with_labels=True, node_size=900, alpha=1.0,
                node_color='cadetblue', font_weight='bold',
                edge_color=['r', 'k', 'r', 'k', 'k'])
        # nx.draw_networkx_edge_labels(graph, pos=pos)
        plt.text(2, 0.1, 'non-linear interaction in red')
        # see here for an example on how to plot edge labels:
        # http://stackoverflow.com/questions/10104700/how-to-set-networkx-edge-labels-offset-to-avoid-label-overlap
        plt.show()
        return fig


if __name__ == "__main__":
    p = undirected()
    print(p.adaa())
    print(p.plot_mute_graph())