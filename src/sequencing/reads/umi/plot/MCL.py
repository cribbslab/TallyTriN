

class mcl(object):

    def __init__(self, ):
        pass

    def example(self, ):
        import markov_clustering as mc
        import networkx as nx
        import random

        # number of nodes to use
        numnodes = 300

        # generate random positions as a dictionary where the key is the node id and the value
        # is a tuple containing 2D coordinates
        positions = {i: (random.random() * 2 - 0, random.random() * 2 +8) for i in range(numnodes)}
        print(positions)
        # use networkx to generate the graph
        network = nx.random_geometric_graph(numnodes, radius=0.3, pos=positions)
        # network = nx.thresholded_random_geometric_graph(numnodes, theta=1.9, radius=0.6, pos=positions)
        # network = nx.random_tree(numnodes)
        # then get the adjacency matrix (in sparse form)
        matrix = nx.to_scipy_sparse_matrix(network)
        print(matrix)

        result = mc.run_mcl(matrix)  # run MCL with default parameters
        clusters = mc.get_clusters(result)  # get clusters

        mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")
        import matplotlib.pyplot as plt
        plt.annotate('UMI number=' + str(numnodes), (-0.4, 1.1))

        plt.show()
        result = mc.run_mcl(matrix, inflation=1.4)
        clusters = mc.get_clusters(result)
        mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")

        # cluster using the optimized cluster inflation value
        result = mc.run_mcl(matrix, inflation=2.1)
        clusters = mc.get_clusters(result)

        mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")

        return


if __name__ == "__main__":
    p = mcl()

    print(p.example())