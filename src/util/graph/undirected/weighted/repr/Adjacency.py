__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.util.sequence.fastq.Read import read as rfastq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to
import numpy as np
import pandas as pd
from src.sequencing.reads.umi.Filter import filter
from src.sequencing.reads.umi.RuleOut import ruleOut as umiro
from src.sequencing.reads.seq.RuleOut import ruleOut as seqro
from src.sequencing.reads.similarity.distance.Hamming import hamming
import networkx as nx
import matplotlib.pyplot as plt
from collections import deque


class adjacency(object):

    def __init__(self, graph):
        self._graph = graph
        self._graph_mapped = None
        self.glen = len(self._graph)

    @property
    def graph(self, ):
        print('the current graph is {}'.format(self._graph))
        return self._graph

    @graph.setter
    def graph(self, value):
        self._graph = value

    @property
    def graph_mapped(self, ):
        if self._graph == None:
            raise 'go set your graph'
        else:
            return self.map(self._graph)

    @property
    def key_mapped(self, ):
        return {[*self._graph.keys()][k]: k for k in range(len([*self._graph.keys()]))}

    def map(self, graph):
        print('->the graph is being mapped')
        print('--->key map: {}'.format(self.key_mapped))
        if isinstance(graph, dict):
            print('--->the graph is a dict')
            g_mapped = {}
            for k, vals in self._graph.items():
                g_mapped[self.key_mapped[k]] = []
                for val in vals:
                    g_mapped[self.key_mapped[k]].append(self.key_mapped[val])
            print('--->the mapped graph: {}'.format(g_mapped))
            return g_mapped

    def dict(self, ):
        return self._graph

    def set(self, ):
        adj_set = {}
        for k, vals in self._graph.items():
            adj_set[k] = set(vals)
        return adj_set

    def list(self, ):
        return [*self._graph.values()]

    def matrix(self, ):
        adj_mat = np.zeros(shape=[self.glen, self.glen])
        for k, vals in self._graph.items():
            for val in vals:
                adj_mat[self.key_mapped[k], self.key_mapped[val]] = 1
        return adj_mat

    def hash(self, ):
        return

    def toEdgeList(self, rr=True):
        edges = []
        for k, vals in self._graph.items():
            for val in vals:
                edges.append((k, val))
        if rr:
            repeat = []
            edge_set = set(edges)
            while edges:
                edge = edges.pop(0)
                if tuple(reversed(edge)) in edges:
                    repeat.append(edge)
            edges = list(edge_set.difference(set(repeat)))
        return edges


if __name__ == "__main__":
    graph_adj_dict = {
        'A': ['B', 'C'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['E', 'F'],
        'E': ['D'],
        'F': ['D'],
    }

    # graph_adj_dict = {
    #     0: [4],
    #     1: [2, 5],
    #     2: [1, 3, 8],
    #     3: [2, 14],
    #     4: [0, 9],
    #     5: [1],
    #     6: [10, 12],
    #     7: [13],
    #     8: [2, 14],
    #     9: [4, 10, 15],
    #     10: [6, 9],
    #     11: [17],
    #     12: [6],
    #     13: [7, 19, 20],
    #     14: [3, 8],
    #     15: [9, 21],
    #     16: [22],
    #     17: [11, 18],
    #     18: [17, 19],
    #     19: [13, 18, 26],
    #     20: [13, 26],
    #     21: [15, 27],
    #     22: [16, 23],
    #     23: [22, 24, 28],
    #     24: [23, 25, 29],
    #     25: [24],
    #     26: [19, 20, 30, 31],
    #     27: [21],
    #     28: [23, 29],
    #     29: [24, 28],
    #     30: [26, 31],
    #     31: [26, 30],
    # }

    # graph = {
    #     0: [(0, 1), (0, 2), (0, 3)],
    #     1: [],
    #     2: [(2, 1)],
    #     3: [(3, 4), (3, 5)],
    #     4: [(4, 3), (4, 5)],
    #     5: [(5, 3), (5, 4), (5, 7)],
    #     6: [(6, 8)],
    #     7: [],
    #     8: [(8, 9)],
    #     9: []
    # }

    p = adjacency(graph_adj_dict)

    p.graph = p.graph_mapped

    print(p.graph)

    print(p.dict())

    print(p.set())

    print(p.list())

    print(p.matrix())

    print(p.toEdgeList(rr=False))
    print(p.toEdgeList(rr=True))