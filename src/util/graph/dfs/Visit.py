from collections import OrderedDict


class visit(object):

    def __init__(self, graph):
        self.graph = graph

    def recursive(self, node):
        visited = []
        def recur(node):
            visited.append(node)
            print(visited)
            for neighbor in self.graph[node]:
                print(neighbor)
                if neighbor not in visited:
                    recur(neighbor)
        recur(node)
        return visited

    def d(self, node, node_val_sorted, node_set_remaining):
        visited = set()
        def recur(node):
            visited.add(node)
            # print(visited)
            for neighbor in self.graph[node]:
                # print(neighbor)
                if neighbor not in visited:
                    if neighbor in node_set_remaining:
                        if node_val_sorted[node] >= 2 * node_val_sorted[neighbor] - 1:
                            recur(neighbor)
        recur(node)
        return visited


if __name__ == "__main__":
    graph = {
        'A': ['B', 'C'],
        'B': ['A', 'D'],
        'C': ['A'],
        'D': ['B', 'E', 'F'],
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

    node_val_sorted = {
        'A': 456,
        'E': 90,
        'D': 72,
        'B': 2,
        'C': 2,
        'F': 1,
    }

    p = visit(graph)
    # print(p.recursive('A'))

    nodes = [*node_val_sorted.keys()]
    print(nodes)
    node_cp = nodes.copy()
    node_set_remaining = set(node_cp)
    cc_sub = []
    while nodes:
        e = nodes.pop(0)
        if e in node_set_remaining:
            seen = p.d(e, node_val_sorted, node_set_remaining)
            cc_sub.append(list(seen))
            node_set_remaining = node_set_remaining - seen
            print('remaining: {}'.format(node_set_remaining))
        else:
            continue
    print(cc_sub)